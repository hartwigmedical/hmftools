package com.hartwig.hmftools.common.genome.refgenome;

public enum RefGenomeVersion
{
    HG19,
    HG37, // included to distinguish from HG19 since has a 'chr' prefix
    HG38;

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";

    public static final String CHR_PREFIX = "chr";

    public static String refGenomeChromosome(final String chromosome, RefGenomeVersion version)
    {
        if((version == HG38 || version == HG37) && !chromosome.contains(CHR_PREFIX))
            return CHR_PREFIX + chromosome;
        else if(version == HG19)
            return stripChromosome(chromosome);
        else
            return chromosome;
    }

    public static String stripChromosome(final String chromosome)
    {
        if(chromosome.startsWith(CHR_PREFIX))
            return chromosome.substring(CHR_PREFIX.length());

        return chromosome;
    }

}
