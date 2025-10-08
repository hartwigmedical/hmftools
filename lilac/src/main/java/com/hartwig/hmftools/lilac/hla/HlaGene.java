package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.GeneClass.MHC_CLASS_1;
import static com.hartwig.hmftools.lilac.GeneClass.MHC_CLASS_2;
import static com.hartwig.hmftools.lilac.GeneClass.PGX;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;

import com.hartwig.hmftools.lilac.GeneClass;

public enum HlaGene
{
    HLA_A(MHC_CLASS_1, false),
    HLA_B(MHC_CLASS_1, false),
    HLA_C(MHC_CLASS_1, false),

    HLA_Y(MHC_CLASS_1, true),
    HLA_H(MHC_CLASS_1, true),

    HLA_DQB1(MHC_CLASS_2, false),
    HLA_DPA1(MHC_CLASS_2, false),
    HLA_DPB1(MHC_CLASS_2, false),
    HLA_DQA1(MHC_CLASS_2, false),
    HLA_DRB1(MHC_CLASS_2, false),
    HLA_DRB3(MHC_CLASS_2, false, false, false),
    HLA_DRB4(MHC_CLASS_2, false, false, false),
    HLA_DRB5(MHC_CLASS_2, false, false, false),

    DPYD(PGX, false),

    NONE(MHC_CLASS_1, false, true); // used for debugging

    private final GeneClass mMhcClass;
    private final boolean mIsPseudo;
    private final boolean mIsDebug;
    private final boolean mHasFrequencies;

    HlaGene(final GeneClass mhcClass, final boolean isPseudo)
    {
        this(mhcClass, isPseudo, false, true);
    }

    HlaGene(final GeneClass mhcClass, final boolean isPseudo, final boolean isDebug)
    {
        this(mhcClass, isPseudo, isDebug, true);
    }

    HlaGene(final GeneClass mhcClass, final boolean isPseudo, final boolean isDebug, final boolean hasFrequencies)
    {
        mMhcClass = mhcClass;
        mIsPseudo = isPseudo;
        mIsDebug = isDebug;
        mHasFrequencies = hasFrequencies;
    }

    public GeneClass mhcClass() { return mMhcClass; }
    public boolean isPseudo() { return mIsPseudo; }
    public boolean isDebug() { return mIsDebug; }
    public boolean hasFrequencies() { return mHasFrequencies; }

    public static HlaGene fromString(final String s)
    {
        if(s.equals("DPYD"))
            return DPYD;

        String geneStr = s.startsWith(HLA_PREFIX) ? s.substring(HLA_PREFIX.length()) : s;
        try
        {
            HlaGene gene = valueOf("HLA_" + geneStr);
            return gene;
        }
        catch(IllegalArgumentException e)
        {
            return null;
        }
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

        if(mMhcClass == PGX)
            return longName();

        return super.toString().substring(HLA_PREFIX.length());
    }
}
