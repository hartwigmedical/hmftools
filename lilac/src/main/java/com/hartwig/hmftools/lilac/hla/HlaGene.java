package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;

public enum HlaGene
{
    HLA_A,
    HLA_B,
    HLA_C,
    HLA_Y,
    HLA_H,
    NONE; // used for debugging

    public static HlaGene fromString(final String s)
    {
        String gene = s.startsWith(HLA_PREFIX) ? s.substring(HLA_PREFIX.length()) : s;
        return valueOf("HLA_" + gene);
    }

    @Override
    public String toString()
    {
        return longName();
    }

    public String longName()
    {
        if(this == NONE)
            return "";

        return super.toString().replace("_", "-");
    }

    public String shortName()
    {
        if(this == NONE)
            return "";

        return super.toString().substring(HLA_PREFIX.length());
    }
}
