package com.hartwig.hmftools.cup.traits;

public enum SampleTraitType
{
    GENDER,
    WGD,
    PURITY,
    PLOIDY,
    SNV_COUNT,
    MS_INDELS_TMB,
    CHORD_HRD;

    public String getAlias()
    {
        switch(this) {
            case GENDER:
                return "is_male";
            case WGD:
                return "whole_genome_duplication";
            case PURITY:
                return "purity";
            case PLOIDY:
                return "ploidy";
            case SNV_COUNT:
                return "snv_count";
            case MS_INDELS_TMB:
                return "indels_per_mb";
            case CHORD_HRD:
                return "chord_hrd";
            default:
                throw new IllegalArgumentException("Alias not implemented for: " + name());
        }
    }
}
