package com.hartwig.hmftools.common.genome.chromosome;

public interface Chromosome {

    String contig();

    boolean isAutosome();

    boolean isAllosome();

    static boolean isAltRegionContig(final String regionContig)
    {
        return !HumanChromosome.contains(regionContig) && !MitochondrialChromosome.contains(regionContig);
    }
}
