package com.hartwig.hmftools.cdr3;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class TestVJGeneStore implements VJGeneStore
{
    private final List<VJGene> mVJGenes;
    private final Multimap<String, VJGene> mAnchorSequenceMap = ArrayListMultimap.create();

    private final Map<VJGeneType, Multimap<String, VJGene>> mGeneTypeAnchorSeqMap = new HashMap<>();
    private final Multimap<VJAnchorReferenceLocation, VJGene> mGeneLocationVJGeneMap = ArrayListMultimap.create();

    @Override
    public List<VJGene> getVJGenes()
    {
        return mVJGenes;
    }

    @Override
    public Set<String> getAnchorSequenceSet()
    {
        return mAnchorSequenceMap.keySet();
    }

    @Override
    public Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType)
    {
        Multimap<String, VJGene> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.keySet() : Collections.emptySet();
    }

    @Override
    public Collection<VJGene> getByAnchorSequence(@NotNull String anchorSeq)
    {
        return mAnchorSequenceMap.get(anchorSeq);
    }

    @Override
    public Collection<VJGene> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq)
    {
        Multimap<String, VJGene> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.get(anchorSeq) : Collections.emptySet();
    }

    @Override
    public Collection<VJGene> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation)
    {
        return mGeneLocationVJGeneMap.get(vjAnchorReferenceLocation);
    }

    @Override
    public Collection<VJAnchorReferenceLocation> getVJAnchorReferenceLocations()
    {
        return mGeneLocationVJGeneMap.keySet();
    }

    public TestVJGeneStore(List<VJGene> vjGenes)
    {
        mVJGenes = vjGenes;

        // from this we find all the anchor sequence locations and fix them
        for (VJGene gene : mVJGenes)
        {
            if (gene.getAnchorLocation() != null)
            {
                mGeneLocationVJGeneMap.put(new VJAnchorReferenceLocation(gene.getType().getVj(), gene.getAnchorLocation()), gene);
            }

            if (!gene.getAnchorSequence().isEmpty())
            {
                mAnchorSequenceMap.put(gene.getAnchorSequence(), gene);
                mGeneTypeAnchorSeqMap.computeIfAbsent(gene.getType(), o -> ArrayListMultimap.create())
                        .put(gene.getAnchorSequence(), gene);
            }
        }
    }
}