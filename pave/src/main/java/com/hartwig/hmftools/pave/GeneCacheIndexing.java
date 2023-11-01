package com.hartwig.hmftools.pave;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;

public class GeneCacheIndexing
{
    public final List<GeneData> ChromosomeGenes;
    public final List<GeneData> CurrentGenes; // in the current vacinity
    public int CurrentPosStrandGeneIndex;
    public int CurrentNegStrandGeneIndex;

    public GeneCacheIndexing(final List<GeneData> genes)
    {
        ChromosomeGenes = genes;
        CurrentGenes = Lists.newArrayList();
        CurrentNegStrandGeneIndex = 0;
        CurrentPosStrandGeneIndex = 0;
    }
}
