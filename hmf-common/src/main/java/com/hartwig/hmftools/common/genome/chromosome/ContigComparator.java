package com.hartwig.hmftools.common.genome.chromosome;

import java.util.Comparator;

public enum ContigComparator implements Comparator<String> {

    INSTANCE;

    @Override
    public int compare(final String contig1, final String contig2) {
        final boolean humanContig1 = HumanChromosome.contains(contig1);
        final boolean humanContig2 = HumanChromosome.contains(contig2);

        if (humanContig1 || humanContig2) {
            if (humanContig1 && humanContig2) {
                return HumanChromosome.fromString(contig1).compareTo(HumanChromosome.fromString(contig2));
            }
            return humanContig1 ? -1 : 1;
        }

        final boolean mitochondrialContig1 = MitochondrialChromosome.contains(contig1);
        final boolean mitochondrialContig2 = MitochondrialChromosome.contains(contig2);

        if (mitochondrialContig1 || mitochondrialContig2) {
            if (mitochondrialContig1 && mitochondrialContig2) {
                return MitochondrialChromosome.fromString(contig1).compareTo(MitochondrialChromosome.fromString(contig2));
            }
            return mitochondrialContig1 ? -1 : 1;
        }

        return contig1.compareTo(contig2);
    }
}
