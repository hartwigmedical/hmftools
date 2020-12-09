package com.hartwig.hmftools.common.genome.refgenome;

public enum RefGenomeVersion
{
    RG_37,
    RG_38,
    RG_19; // included to distinguish from HG37 since has has the 'chr' prefix

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";

    public static final String CHR_PREFIX = "chr";

    public RefGenomeVersion defaultVersion() { return RG_37; }

    public static boolean is37(final RefGenomeVersion version) { return version == RG_37 || version == RG_19; }

    public static RefGenomeVersion from(final String version)
    {
        if(version.equals(RG_37.toString()) || version.equals("37") || version.equals("HG37"))
            return RG_37;

        if(version.equals(RG_38.toString()) || version.equals("38") || version.equals("HG38"))
            return RG_38;

        if(version.equals(RG_19.toString()) || version.equals("19") || version.equals("HG19"))
            return RG_19;

        return RG_37;
    }

    public static String refGenomeChromosome(final String chromosome, RefGenomeVersion version)
    {
        if((version == RG_38 || version == RG_19) && !chromosome.contains(CHR_PREFIX))
            return CHR_PREFIX + chromosome;
        else if(version == RG_37)
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
