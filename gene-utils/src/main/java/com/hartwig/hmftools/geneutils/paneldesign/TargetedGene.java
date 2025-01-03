package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;

public class TargetedGene
{
    private final GeneData mGeneData;
    private final TranscriptData mTranscriptData;

    // regions, sorted by position
    private final List<TargetedGeneRegion> mRegions = new ArrayList<>();

    public GeneData getGeneData()
    {
        return mGeneData;
    }

    public TranscriptData getTranscriptData()
    {
        return mTranscriptData;
    }

    public List<TargetedGeneRegion> getRegions()
    {
        return mRegions;
    }

    public TargetedGene(final GeneData geneData, final TranscriptData transcriptData)
    {
        mGeneData = geneData;
        mTranscriptData = transcriptData;
    }

    public void addRegion(TargetedGeneRegion.Type type, int positionStart, int positionEnd)
    {
        mRegions.add(new TargetedGeneRegion(this, type, new BaseRegion(positionStart, positionEnd)));
    }
}
