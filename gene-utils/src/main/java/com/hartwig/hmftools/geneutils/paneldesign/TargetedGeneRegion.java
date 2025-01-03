package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

public class TargetedGeneRegion
{
    public enum Type
    {
        CODING, UTR, UP_STREAM, DOWN_STREAM, INTRONIC_LONG, INTRONIC_SHORT
    }

    private final TargetedGene mGene;
    private final Type mType;
    private final BaseRegion mBaseRegion;

    private ProbeCandidate mSelectedProbe = null;

    private final List<ProbeCandidate> mProbeCandidates = new ArrayList<>();

    public TargetedGeneRegion(final TargetedGene gene, Type type, BaseRegion region)
    {
        mGene = gene;
        mType = type;
        mBaseRegion = region;
    }

    public TargetedGene getGene()
    {
        return mGene;
    }

    public Type getType()
    {
        return mType;
    }

    public String getChromosome()
    {
        return mGene.getGeneData().Chromosome;
    }

    public int getStart()
    {
        return mBaseRegion.start();
    }

    public int getEnd()
    {
        return mBaseRegion.end();
    }

    public ProbeCandidate getSelectedProbe()
    {
        return mSelectedProbe;
    }

    public void setSelectedProbe(final ProbeCandidate probe)
    {
        mSelectedProbe = probe;
    }

    public List<ProbeCandidate> getProbeCandidates()
    {
        return mProbeCandidates;
    }

    public boolean useWholeRegion()
    {
        return mType == Type.CODING || mType == Type.UTR;
    }
}
