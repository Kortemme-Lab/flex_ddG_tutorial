#!/usr/bin/env python3
"""
Generate a copy of ddG-backrub.xml that additionally reports per-chain energies.

Why
---
The `unbound_wt` / `unbound_mut` totals in ddG.db3 are whole-pose scores. InterfaceDdGMover
builds the unbound state by rigid-body translating the moving chain(s) 1000 A apart, without
repacking or re-minimizing, and every chain stays in the pose
(protocols/features/InterfaceDdGMover.cc::unbind). So

    unbound_X_total = sum over chains of intra(chain)        for X in {wt, mut}

and `unbound_mut - unbound_wt` is the summed intramolecular ddG of *every* chain, including
the ones that were never mutated. Those only differ between the wild type and mutant branches
because the whole pose is re-minimized independently in each, so they contribute artifact.

What this adds
--------------
A `TotalEnergyMetric` per chain, reported into ddG.db3 by the same ReportToDB mover that
InterfaceDdGMover already applies to all four states. This is *reporting only* -- it adds no
sampling, consumes no random numbers, and leaves ddG_bind bit-identical (verified against the
unmodified protocol with -constant_seed).

Reading the result
------------------
On the **unbound** states the chains are 1000 A apart, so there are no cross-chain pair
energies and the reported value is exactly that chain's intramolecular energy. That is the
number to use:

    per-chain ddG(chain X) = intra_mut(X) - intra_wt(X)

On the **bound** states the value is intra(chain) + roughly half the interface energy, because
Rosetta's residue_total_energies splits each two-body term between its two residues. Those rows
are useful as a cross-check (they sum to the pose total) but should not be read as per-chain
stability numbers.

Caveat this does NOT fix
------------------------
The unbound state is still never relaxed. This gives you the intramolecular strain difference
*in the bound backbone conformation*, not a folding ddG of the free monomer. For that you want
a dedicated monomer protocol (e.g. cartesian_ddg) run on the isolated chain.
"""

import os


def chains_in_pdb(pdb_path):
    """Chain IDs in the order they first appear in the PDB."""
    chains = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21]
                if chain not in chains:
                    chains.append(chain)
    return chains


def write_per_chain_protocol(base_xml_path, chains, out_path, scorefxn='fa_talaris2014'):
    """Write a copy of base_xml_path with a per-chain TotalEnergyMetric for each chain.

    Each chain gets its own database table (chain_<X>_energy) with a single `total_energy`
    column. Separate tables rather than one table with prefixed columns is deliberate: as of
    Rosetta 2022.45, SimpleMetricFeatures builds schema column names as `custom_type + name`
    but builds the INSERT as `custom_type + "_" + name`, so any non-empty custom_type produces
    "table simple_metrics has no column named ..." at report time.
    """
    if not chains:
        raise ValueError('No chains given')

    # Rosetta lowercases table names, so chains differing only in case would collide.
    lowered = [c.lower() for c in chains]
    if len(set(lowered)) != len(lowered):
        raise ValueError('Chain IDs differ only by case, which collides in table names: %s' % chains)

    xml = open(base_xml_path).read()

    def substitute(text, anchor, addition):
        if anchor not in text:
            raise ValueError('Could not find anchor in %s:\n  %s' % (base_xml_path, anchor))
        return text.replace(anchor, anchor + addition, 1)

    selectors = ''.join(
        '\n    <Chain name="chain_%s" chains="%s"/>' % (c, c) for c in chains)
    xml = substitute(
        xml, '    <StoredResidueSubset name="restore_neighbor_shell" subset_name="neighbor_shell"/>',
        selectors)

    metrics = '\n'.join(
        '    <TotalEnergyMetric name="chain_%s_energy" residue_selector="chain_%s" scorefxn="%s"/>'
        % (c, c, scorefxn) for c in chains)
    anchor = '  <FILTERS>'
    if anchor not in xml:
        raise ValueError('Could not find FILTERS block in %s' % base_xml_path)
    xml = xml.replace(anchor, '  <SIMPLE_METRICS>\n%s\n  </SIMPLE_METRICS>\n\n%s' % (metrics, anchor), 1)

    reporters = ''.join(
        '\n      <SimpleMetricFeatures metrics="chain_%s_energy" table_name="chain_%s_energy"/>'
        % (c, c) for c in chains)
    xml = substitute(
        xml, '      <StructureScoresFeatures scorefxn="%s"/>' % scorefxn, reporters)

    with open(out_path, 'w') as f:
        f.write(xml)
    return out_path


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 4:
        sys.exit('usage: per_chain_protocol.py <base.xml> <input.pdb> <out.xml>')
    base_xml, pdb, out = sys.argv[1:4]
    found = chains_in_pdb(pdb)
    write_per_chain_protocol(base_xml, found, out)
    print('Wrote %s with per-chain metrics for chains: %s' % (out, ', '.join(found)))
