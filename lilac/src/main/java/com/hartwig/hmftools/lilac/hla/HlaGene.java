package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_1;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_2;

import com.hartwig.hmftools.lilac.MhcClass_;

public enum HlaGene
{
    HLA_A(CLASS_1, false),
    HLA_B(CLASS_1, false),
    HLA_C(CLASS_1, false),

    HLA_Y(CLASS_1, true),
    HLA_H(CLASS_1, true),

    HLA_DQB1(CLASS_2, false),

    NONE(CLASS_1, false); // used for debugging

    private final MhcClass_ mMhcClass;
    private final boolean mIsPseudo;

    HlaGene(final MhcClass_ mhcClass, boolean isPseudo)
    {
        mMhcClass = mhcClass;
        mIsPseudo = isPseudo;
    }

    public MhcClass_ mhcClass()
    {
        return mMhcClass;
    }

    public boolean isPseudo()
    {
        return mIsPseudo;
    }

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
