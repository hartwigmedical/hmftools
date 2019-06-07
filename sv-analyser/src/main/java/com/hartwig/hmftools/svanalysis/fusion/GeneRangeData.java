package com.hartwig.hmftools.svanalysis.fusion;

import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;

public class GeneRangeData
{
    public static int GENE_PHASING_REGION_5P_UTR = 0;
    public static int GENE_PHASING_REGION_CODING_0 = 1;
    public static int GENE_PHASING_REGION_CODING_1 = 2;
    public static int GENE_PHASING_REGION_CODING_2 = 3;
    public static int GENE_PHASING_REGION_PROMOTOR = 4;
    public static int GENE_PHASING_REGION_MAX = 5;

    public final EnsemblGeneData GeneData;
    public final String Arm;
    private int[] mFivePrimePhaseCounts;
    private int[] mThreePrimePhaseCounts;

    public GeneRangeData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;
        mFivePrimePhaseCounts = new int[GENE_PHASING_REGION_MAX];
        mThreePrimePhaseCounts = new int[GENE_PHASING_REGION_MAX];

        Arm = SvUtilities.getChromosomalArm(geneData.Chromosome, geneData.GeneStart);
    }

    public int[] getPhaseCounts(boolean useFive) { return useFive ? mFivePrimePhaseCounts : mThreePrimePhaseCounts; }

    public boolean hasCodingTranscripts(boolean useFive)
    {
        return hasCodingTranscripts(getPhaseCounts(useFive));
    }

    public boolean hasCodingTranscripts(final int[] phaseCounts)
    {
        return phaseCounts[GENE_PHASING_REGION_CODING_0] > 0
            || phaseCounts[GENE_PHASING_REGION_CODING_1] > 0
            || phaseCounts[GENE_PHASING_REGION_CODING_2] > 0;
    }

}
