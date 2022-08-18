package com.hartwig.hmftools.cider;

import java.util.ArrayList;
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
    private final Multimap<String, VJAnchorTemplate> mAnchorSequenceMap = ArrayListMultimap.create();

    private final Map<VJGeneType, Multimap<String, VJAnchorTemplate>> mGeneTypeAnchorSeqMap = new HashMap<>();
    private final Multimap<VJAnchorReferenceLocation, VJAnchorTemplate> mGeneLocationVJGeneMap = ArrayListMultimap.create();

    private final List<IgTcrConstantRegion> mIgTcrConstantRegions = new ArrayList<>();

    @Override
    public Set<String> getAnchorSequenceSet()
    {
        return mAnchorSequenceMap.keySet();
    }

    @Override
    public Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType)
    {
        Multimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.keySet() : Collections.emptySet();
    }

    @Override
    public Collection<VJAnchorTemplate> getByAnchorSequence(@NotNull String anchorSeq)
    {
        return mAnchorSequenceMap.get(anchorSeq);
    }

    @Override
    public Collection<VJAnchorTemplate> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq)
    {
        Multimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.get(anchorSeq) : Collections.emptySet();
    }

    @Override
    public Collection<VJAnchorTemplate> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation)
    {
        return mGeneLocationVJGeneMap.get(vjAnchorReferenceLocation);
    }

    @Override
    public Collection<VJAnchorReferenceLocation> getVJAnchorReferenceLocations()
    {
        return mGeneLocationVJGeneMap.keySet();
    }

    @Override
    public Collection<IgTcrConstantRegion> getIgConstantRegions() { return mIgTcrConstantRegions; }

    public TestVJGeneStore(List<VJAnchorTemplate> vjAnchorTemplates)
    {
        // from this we find all the anchor sequence locations and fix them
        for (VJAnchorTemplate gene : vjAnchorTemplates)
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