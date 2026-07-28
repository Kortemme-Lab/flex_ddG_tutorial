#!/usr/bin/env python3
"""
Recover per-chain intramolecular energies from an EXISTING flex ddG run.

No re-sampling is required. The backrub trajectory is the expensive part and you already ran
it: struct.db3 holds the coordinates of every backrub, wild type minimized and mutant minimized
pose. This script reads those poses back into Rosetta, isolates one chain at a time, rescores,
and reports per-chain ddG. Verified to reproduce the in-protocol per-chain metric (see
per_chain_protocol.py) to nine decimal places.

Why the chains have to be physically isolated rather than just selected: struct.db3 stores the
*bound* poses. A Chain residue selector scoped over a bound pose still picks up cross-chain pair
energies, because Rosetta's residue_total_energies splits each two-body term between its two
residues, so roughly half the interface energy leaks into each chain. Deleting the other chains
reproduces the separated unbound state exactly, since intra-chain energy is invariant under the
rigid-body translation InterfaceDdGMover uses to unbind.

What this does and does not give you
------------------------------------
It gives you the intramolecular energy difference between the mutant and wild type in the
*bound* backbone conformation, per chain. That is a strain term. It is NOT a folding ddG of the
free monomer: the unbound state is never relaxed here, in this script or in flex ddG itself. For
a true monomer stability ddG use a dedicated protocol such as cartesian_ddg on the isolated chain.

Built-in control: a chain you did not mutate should come out at 0 within its SEM. If it does not,
nstruct is too low to average out the whole-pose minimization noise, and the mutated chain's
number is not trustworthy either.

Usage
-----
    python3 reprocess_per_chain.py <output_folder> [--stride N] [--chains A,B] [--csv out.csv]

--stride must match the backrub_trajectory_stride you ran with (default 5, as in the tutorial).
It only affects the checkpoint labels, not any energy.
"""

import argparse
import glob
import os
import re
import sqlite3
import subprocess
import sys

import numpy as np
import pandas as pd

rosetta_scripts_path = os.path.expanduser('~/rosetta/source/bin/rosetta_scripts')

# Must match the score function the original run used.
rosetta_flags = [
    '-restore_talaris_behavior',
    '-in:file:fullatom',
    '-out:nooutput',
]

struct_db3_name = 'struct.db3'
per_chain_db3_name = 'per_chain.db3'

# flex ddG writes three poses per checkpoint, in this order.
pose_order = ['backrub', 'wt', 'mut']


def chains_in_struct_db3(struct_db3):
    conn = sqlite3.connect(struct_db3)
    try:
        chains = [row[0] for row in conn.execute(
            'SELECT DISTINCT chain_id FROM residue_pdb_identification ORDER BY chain_id')]
    finally:
        conn.close()
    return [c for c in chains if c and c.strip()]


def write_rescore_protocol(chains, out_path, scorefxn='fa_talaris2014'):
    """Emit a RosettaScripts protocol that isolates and rescores each chain in turn."""
    lowered = [c.lower() for c in chains]
    if len(set(lowered)) != len(lowered):
        raise ValueError('Chain IDs differ only by case, which collides in batch names: %s' % chains)

    selectors = '\n'.join(
        '    <Chain name="chain_%s" chains="%s"/>\n'
        '    <Not name="not_chain_%s" selector="chain_%s"/>' % (c, c, c, c) for c in chains)

    movers = '\n'.join(
        '    <DeleteRegionMover name="isolate_chain_%s" residue_selector="not_chain_%s"/>\n'
        '    <ReportToDB name="chain_%s_report" batch_description="per_chain" database_name="%s">\n'
        '      <ScoreTypeFeatures/>\n'
        '      <ScoreFunctionFeatures scorefxn="%s"/>\n'
        '      <StructureScoresFeatures scorefxn="%s"/>\n'
        '    </ReportToDB>' % (c, c, c, per_chain_db3_name, scorefxn, scorefxn) for c in chains)

    steps = ['    <Add mover_name="save_full"/>']
    for i, c in enumerate(chains):
        if i > 0:
            steps.append('    <Add mover_name="restore_full"/>')
        steps.append('    <Add mover_name="isolate_chain_%s"/>' % c)
        steps.append('    <Add mover_name="chain_%s_report"/>' % c)

    xml = '''<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="%s" weights="talaris2014"/>
  </SCOREFXNS>

  <RESIDUE_SELECTORS>
%s
  </RESIDUE_SELECTORS>

  <MOVERS>
    <SavePoseMover name="save_full" reference_name="full_pose" restore_pose="0"/>
    <SavePoseMover name="restore_full" reference_name="full_pose" restore_pose="1"/>
%s
  </MOVERS>

  <PROTOCOLS>
%s
  </PROTOCOLS>
  <OUTPUT />
</ROSETTASCRIPTS>
''' % (scorefxn, selectors, movers, '\n'.join(steps))

    with open(out_path, 'w') as f:
        f.write(xml)
    return out_path


def rescore_one(struct_db3, protocol_path):
    """Run the isolate-and-rescore protocol on one struct.db3, writing per_chain.db3 beside it."""
    working_dir = os.path.dirname(os.path.abspath(struct_db3))
    out_db3 = os.path.join(working_dir, per_chain_db3_name)
    if os.path.isfile(out_db3):
        os.remove(out_db3)

    args = [
        os.path.abspath(rosetta_scripts_path),
        '-inout:dbms:database_name', struct_db3_name,
        '-in:use_database',
        '-parser:protocol', os.path.abspath(protocol_path),
    ] + rosetta_flags

    log_path = os.path.join(working_dir, 'per_chain_rescore.log')
    with open(log_path, 'w') as log:
        proc = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT, cwd=working_dir)
        returncode = proc.wait()

    if returncode != 0 or not os.path.isfile(out_db3):
        raise RuntimeError('Rescoring failed for %s -- see %s' % (struct_db3, log_path))
    return out_db3


def read_per_chain_db3(per_chain_db3, struct_number, case_name, stride):
    conn = sqlite3.connect(per_chain_db3)
    df = pd.read_sql_query('''
    SELECT batches.name AS batch, structures.tag AS tag, structure_scores.score_value AS energy
    FROM structure_scores
    INNER JOIN structures ON structures.struct_id=structure_scores.struct_id
    INNER JOIN batches ON batches.batch_id=structure_scores.batch_id
    INNER JOIN score_types ON score_types.batch_id=structure_scores.batch_id
                          AND score_types.score_type_id=structure_scores.score_type_id
    WHERE score_types.score_type_name="total_score"
    ''', conn)
    conn.close()

    # batch name is "chain_<X>_report"; tag is "<original struct.db3 struct_id>_0001"
    df['chain'] = df['batch'].apply(lambda b: b[len('chain_'):-len('_report')])
    original_id = df['tag'].apply(lambda t: int(re.match(r'(\d+)', t).group(1)))
    df['pose'] = original_id.apply(lambda i: pose_order[(i - 1) % len(pose_order)])
    df['backrub_steps'] = original_id.apply(lambda i: stride * (((i - 1) // len(pose_order)) + 1))
    df['struct_num'] = struct_number
    df['case_name'] = case_name
    return df[['case_name', 'struct_num', 'backrub_steps', 'chain', 'pose', 'energy']]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output_folder', help='flex ddG output folder (e.g. "output")')
    parser.add_argument('--stride', type=int, default=5,
                        help='backrub_trajectory_stride used in the original run (default 5)')
    parser.add_argument('--chains', default=None,
                        help='comma-separated chains (default: auto-detect from struct.db3)')
    parser.add_argument('--csv', default=None, help='write the full per-structure table here')
    parser.add_argument('--reuse', action='store_true',
                        help='skip Rosetta where per_chain.db3 already exists')
    args = parser.parse_args()

    if not os.path.isfile(rosetta_scripts_path):
        sys.exit('ERROR: set rosetta_scripts_path to your compiled rosetta_scripts binary')

    struct_db3s = sorted(glob.glob(os.path.join(args.output_folder, '*', '*', struct_db3_name)))
    if not struct_db3s:
        sys.exit('ERROR: no %s found under %s' % (struct_db3_name, args.output_folder))
    print('Found %d %s files' % (len(struct_db3s), struct_db3_name))

    chains = args.chains.split(',') if args.chains else chains_in_struct_db3(struct_db3s[0])
    print('Chains: %s' % ', '.join(chains))

    protocol_path = os.path.join(args.output_folder, 'per_chain_rescore.generated.xml')
    write_rescore_protocol(chains, protocol_path)
    print('Wrote protocol %s\n' % protocol_path)

    frames = []
    for i, struct_db3 in enumerate(struct_db3s, start=1):
        struct_dir = os.path.dirname(struct_db3)
        case_name = os.path.basename(os.path.dirname(struct_dir))
        struct_number = os.path.basename(struct_dir)
        out_db3 = os.path.join(struct_dir, per_chain_db3_name)

        if args.reuse and os.path.isfile(out_db3):
            print('  [%d/%d] %s (reusing)' % (i, len(struct_db3s), struct_dir))
        else:
            print('  [%d/%d] %s' % (i, len(struct_db3s), struct_dir))
            out_db3 = rescore_one(struct_db3, protocol_path)

        frames.append(read_per_chain_db3(out_db3, struct_number, case_name, args.stride))

    per_structure = pd.concat(frames)
    wide = per_structure[per_structure['pose'].isin(['wt', 'mut'])].pivot_table(
        index=['case_name', 'chain', 'backrub_steps', 'struct_num'],
        columns='pose', values='energy').reset_index()
    wide['ddG'] = wide['mut'] - wide['wt']

    if args.csv:
        wide.to_csv(args.csv, index=False)
        print('\nWrote %s' % args.csv)

    summary = wide.groupby(['case_name', 'chain', 'backrub_steps']).agg(
        nstruct=('ddG', 'size'),
        wt_intra=('wt', 'mean'),
        mut_intra=('mut', 'mean'),
        ddG=('ddG', 'mean'),
        ddG_sd=('ddG', 'std'),
    ).reset_index()
    summary['ddG_sem'] = summary['ddG_sd'] / np.sqrt(summary['nstruct'])

    print('\n=== per-chain intramolecular ddG ===')
    print(summary.round(4).to_string(index=False))
    print('\nA chain you did NOT mutate should read ~0 within ddG_sem.')
    print('This is bound-conformation strain, not a folding ddG (see the module docstring).')


if __name__ == '__main__':
    main()
