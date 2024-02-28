package com.hartwig.hmftools.cup.prep;

public enum ItemType
{
    SNV96,
    GEN_POS,
    SIGNATURE,

    TUMOR_MUTATIONAL_BURDEN,
    SAMPLE_TRAIT,
    FUSION,
    DRIVER,
    VIRUS,
    SV_COUNT,

    EXPRESSION,
    ALT_SJ;

    public String getAlias()
    {
        switch(this) {
            case SNV96:
                return "snv96";
            case GEN_POS:
                return "gen_pos";
            case SIGNATURE:
                return "sig";

            case TUMOR_MUTATIONAL_BURDEN:
                return "event.tmb";
            case SAMPLE_TRAIT:
                return "event.trait";
            case FUSION:
                return "event.fusion";
            case DRIVER:
                return "event.driver";
            case VIRUS:
                return "event.virus";
            case SV_COUNT:
                return "event.sv";

            case EXPRESSION:
                return "gene_exp";
            case ALT_SJ:
                return "alt_sj";

            default:
                throw new IllegalArgumentException("Alias not implemented for: " + name());
        }
    }
}
