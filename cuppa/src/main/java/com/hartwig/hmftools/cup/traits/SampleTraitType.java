package com.hartwig.hmftools.cup.traits;

public enum SampleTraitType
{
    GENDER("is_male"),
    WGD("whole_genome_duplication"),
    PURITY("purity"),
    PLOIDY("ploidy"),
    SNV_COUNT("snv_count"),
    MS_INDELS_TMB("indels_per_mb"),
    CHORD_HRD("chord_hrd");

    private final String mAlias;

    SampleTraitType(String alias)
    {
        mAlias = alias;
    }

    public String getAlias()
    {
        return mAlias;
    }
}
