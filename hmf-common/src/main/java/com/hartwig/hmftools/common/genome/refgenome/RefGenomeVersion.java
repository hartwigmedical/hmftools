package com.hartwig.hmftools.common.genome.refgenome;

public enum RefGenomeVersion
{
    HG37,
    HG38;

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";
    public static final int REF_GENOME_HG38 = 38;
    public static final int REF_GENOME_HG37 = 37;

    public static final String CHR_PREFIX = "chr";

    public static String refGenomeChromosome(final String chromosome, RefGenomeVersion version)
    {
        if(version == HG38 && !chromosome.contains(CHR_PREFIX))
            return CHR_PREFIX + chromosome;
        else if(version == HG37)
            return stripChromosome(chromosome);
        else
            return chromosome;
    }

    public static String stripChromosome(final String chromosome)
    {
        return chromosome.startsWith(CHR_PREFIX) ? chromosome.substring(CHR_PREFIX.length()) : chromosome;
    }

}
