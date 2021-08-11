package com.hartwig.hmftools.neo.cohort;

import com.hartwig.hmftools.common.neo.NeoEpitopeType;

public class NeoEpitopeData
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final String GeneName;
    public final double TpmCancer;
    public final double TpmCohort;
    public final int RnaNovelFragments;
    public final int[] RnaBaseDepth;

    public NeoEpitopeData(
            final int id, final NeoEpitopeType variantType, final String variantInfo, final String geneName,
            double tpmCancer, double tpmCohort, int rnaNovelFragments, final int[] rnaBaseDepth)
    {
        Id = id;
        VariantType = variantType;
        VariantInfo = variantInfo;
        GeneName = geneName;
        TpmCancer = tpmCancer;
        TpmCohort = tpmCohort;
        RnaNovelFragments = rnaNovelFragments;
        RnaBaseDepth = rnaBaseDepth;
    }
}
