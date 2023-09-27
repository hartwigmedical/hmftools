package com.hartwig.hmftools.sieve.annotate;

public class RefGenomeRegionAnnotations
{
    public static final String TSV_HEADER = "A Prop\tT Prop\tG Prop\tC Prop";
    private final double mAProp;
    private final double mTProp;
    private final double mGProp;
    private final double mCProp;

    public RefGenomeRegionAnnotations(final double aProp, final double tProp, final double gProp, final double cProp)
    {
        mAProp = aProp;
        mTProp = tProp;
        mGProp = gProp;
        mCProp = cProp;
    }

    public String getTSVFragment()
    {
        return String.valueOf(mAProp) + '\t' + mTProp + '\t' + mGProp + '\t' + mCProp;
    }
}
