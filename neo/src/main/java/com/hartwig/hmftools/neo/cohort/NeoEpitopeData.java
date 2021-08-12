package com.hartwig.hmftools.neo.cohort;

import com.hartwig.hmftools.common.neo.NeoEpitopeType;

public class NeoEpitopeData
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final String GeneId;
    public final String GeneName;
    public final String AminoAcids;
    public final double TpmCancer;
    public final double TpmCohort;
    public final int RnaNovelFragments;
    public final int[] RnaBaseDepth;

    public double GeneExpression;

    public NeoEpitopeData(
            final int id, final NeoEpitopeType variantType, final String variantInfo, final String geneId, final String geneName,
            final String aminoAcids, double tpmCancer, double tpmCohort, int rnaNovelFragments, final int[] rnaBaseDepth)
    {
        Id = id;
        VariantType = variantType;
        VariantInfo = variantInfo;
        GeneId = geneId;
        GeneName = geneName;
        AminoAcids = aminoAcids;
        TpmCancer = tpmCancer;
        TpmCohort = tpmCohort;
        RnaNovelFragments = rnaNovelFragments;
        RnaBaseDepth = rnaBaseDepth;
        GeneExpression = 0;
    }
}
