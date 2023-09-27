package com.hartwig.hmftools.sieve.annotate;

public class RefGenomeRegionAnnotations
{
    public static final String TSV_HEADER = "IsCentromic\tA Prop\tT Prop\tG Prop\tC Prop";

    private final boolean mIsCentromic;
    private final double mAProp;
    private final double mTProp;
    private final double mGProp;
    private final double mCProp;

    public RefGenomeRegionAnnotations(
            final boolean isCentromic,
            final double aProp,
            final double tProp,
            final double gProp,
            final double cProp)
    {
        mIsCentromic = isCentromic;
        mAProp = aProp;
        mTProp = tProp;
        mGProp = gProp;
        mCProp = cProp;
    }

    public String getTSVFragment()
    {
        return String.valueOf(mIsCentromic) + '\t' + mAProp + '\t' + mTProp + '\t' + mGProp + '\t' + mCProp;
    }
}
