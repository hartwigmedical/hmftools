package com.hartwig.hmftools.linx.fusion;

import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;

public class GeneRangeData
{
    public static int GENE_PHASING_REGION_5P_UTR = 0;
    public static int GENE_PHASING_REGION_CODING_0 = 1;
    public static int GENE_PHASING_REGION_CODING_1 = 2;
    public static int GENE_PHASING_REGION_CODING_2 = 3;
    public static int GENE_PHASING_REGION_PROMOTOR = 4;
    public static int GENE_PHASING_REGION_MAX = 5;

    private final EnsemblGeneData GeneData;
    private int[] mPhasingCounts;

    public GeneRangeData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;
        mPhasingCounts = new int[GENE_PHASING_REGION_MAX];
    }

    public final int[] getPhasingCounts() { return mPhasingCounts; }
    public void setPhasingCounts(final int[] counts)
    {
        mPhasingCounts = counts;
    }

}
