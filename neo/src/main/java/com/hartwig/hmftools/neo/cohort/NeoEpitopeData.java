package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;

import java.util.List;

import com.hartwig.hmftools.common.neo.NeoEpitopeType;

public class NeoEpitopeData
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final String GeneId;
    public final String GeneName;
    public final String AminoAcids;
    public final List<String>[] Transcripts;
    public final double TpmCancer;
    public final double TpmCohort;
    public double[] TransExpression;
    public int RnaNovelFragments;
    public final int[] RnaBaseDepth;

    public NeoEpitopeData(
            final int id, final NeoEpitopeType variantType, final String variantInfo, final String geneId, final String geneName,
            final String aminoAcids, final List<String> transNamesUp, final List<String> transNamesDown, double tpmCancer, double tpmCohort)
    {
        Id = id;
        VariantType = variantType;
        VariantInfo = variantInfo;
        GeneId = geneId;
        GeneName = geneName;
        AminoAcids = aminoAcids;
        Transcripts = new List[FS_PAIR];
        Transcripts[FS_UP] = transNamesUp;
        Transcripts[FS_DOWN] = transNamesDown;
        TpmCancer = tpmCancer;
        TpmCohort = tpmCohort;
        RnaNovelFragments = 0;
        RnaBaseDepth = new int[] {0, 0};
        TransExpression = new double[] {-1, -1}; // indicating not set
    }
}
