package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_1;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_2;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.MhcClass_;

public enum HlaGene_
{
    HLA_A(CLASS_1, false),
    HLA_B(CLASS_1, false),
    HLA_C(CLASS_1, false),

    HLA_Y(CLASS_1, true),
    HLA_H(CLASS_1, true),

    HLA_DQB1(CLASS_2, false),
    HLA_DPA1(CLASS_2, false),
    HLA_DPB1(CLASS_2, false),
    HLA_DQA1(CLASS_2, false),

    HLA_DRB1(CLASS_2, false),
    HLA_DRB3(CLASS_2, false, false, false),
    HLA_DRB4(CLASS_2, false, false, false),
    HLA_DRB5(CLASS_2, false, false, false),

    NONE(CLASS_1, false, true); // used for debugging

    private final MhcClass_ mMhcClass;
    private final boolean mIsPseudo;
    private final boolean mIsDebug;
    private final boolean mHasFrequencies;

    HlaGene_(final MhcClass_ mhcClass, boolean isPseudo)
    {
        this(mhcClass, isPseudo, false, true);
    }

    HlaGene_(final MhcClass_ mhcClass, boolean isPseudo, boolean isDebug)
    {
        this(mhcClass, isPseudo, isDebug, true);
    }

    HlaGene_(final MhcClass_ mhcClass, boolean isPseudo, boolean isDebug, boolean hasFrequencies)
    {
        mMhcClass = mhcClass;
        mIsPseudo = isPseudo;
        mIsDebug = isDebug;
        mHasFrequencies = hasFrequencies;
    }

    public MhcClass_ mhcClass()
    {
        return mMhcClass;
    }

    public boolean isPseudo()
    {
        return mIsPseudo;
    }

    public boolean isDebug()
    {
        return mIsDebug;
    }

    public boolean hasFrequencies() { return mHasFrequencies; }

    public static HlaGene_ fromString(final String s)
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
