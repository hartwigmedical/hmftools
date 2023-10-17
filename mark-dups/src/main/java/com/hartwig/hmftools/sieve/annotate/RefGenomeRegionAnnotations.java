package com.hartwig.hmftools.sieve.annotate;

public class RefGenomeRegionAnnotations
{
    public static final String TSV_HEADER = "DistToCentromere\tDistToTelomere\tDistToMasked\tGenes\tA Prop\tT Prop\tG Prop\tC Prop";

    private final int mDistToCentromere;
    private final int mDistToTelomere;
    private final int mDistToMasked;
    private final String mGenesStr;
    private final double mAProp;
    private final double mTProp;
    private final double mGProp;
    private final double mCProp;

    public RefGenomeRegionAnnotations(
            final int distToCentromere,
            final int distToTelomere,
            final int distToMasked,
            final String genesStr,
            final double aProp,
            final double tProp,
            final double gProp,
            final double cProp)
    {
        mDistToCentromere = distToCentromere;
        mDistToTelomere = distToTelomere;
        mDistToMasked = distToMasked;
        mGenesStr = genesStr;
        mAProp = aProp;
        mTProp = tProp;
        mGProp = gProp;
        mCProp = cProp;
    }

    public String getTSVFragment()
    {
        return String.valueOf(mDistToCentromere) + '\t' + mDistToTelomere + '\t' + mDistToMasked + '\t' + mGenesStr + '\t' + mAProp + '\t'
                + mTProp + '\t' + mGProp + '\t' + mCProp;
    }
}
