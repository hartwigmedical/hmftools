package com.hartwig.hmftools.cup.prep;

public enum ItemType
{
    SIGNATURE("sig"),
    SNV96("snv96"),
    GEN_POS("gen_pos"),

    TUMOR_MUTATIONAL_BURDEN("event.tmb"),

    SAMPLE_TRAIT("event.trait"),
    SV_COUNT("event.sv"),
    VIRUS("event.virus"),
    FUSION("event.fusion"),
    DRIVER("event.driver"),

    EXPRESSION("gene_exp"),
    ALT_SJ("alt_sj");

    private final String mAlias;

    ItemType(String alias)
    {
        mAlias = alias;
    }

    public String getAlias()
    {
        return mAlias;
    }

    public static ItemType fromAlias(String alias)
    {
        switch(alias) {
            case "sig": return SIGNATURE;
            case "snv96": return SNV96;
            case "gen_pos": return GEN_POS;
            case "event.tmb": return TUMOR_MUTATIONAL_BURDEN;
            case "event.trait": return SAMPLE_TRAIT;
            case "event.sv": return SV_COUNT;
            case "event.virus": return VIRUS;
            case "event.fusion": return FUSION;
            case "event.driver": return DRIVER;
            case "gene_exp": return EXPRESSION;
            case "alt_sj": return ALT_SJ;
            default:
                throw new IllegalArgumentException("Invalid alias: " + alias);
        }
    }
}
