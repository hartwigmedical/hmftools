package com.hartwig.hmftools.cider;

import java.util.Collections;
import java.util.List;

import org.eclipse.collections.api.collection.ImmutableCollection;
import org.eclipse.collections.api.factory.Lists;
import org.eclipse.collections.api.factory.Sets;
import org.eclipse.collections.api.list.MutableList;
import org.eclipse.collections.api.map.MutableMap;
import org.eclipse.collections.api.multimap.Multimap;
import org.eclipse.collections.api.multimap.MutableMultimap;
import org.eclipse.collections.api.set.SetIterable;
import org.eclipse.collections.impl.list.mutable.FastList;
import org.eclipse.collections.impl.map.mutable.UnifiedMap;
import org.eclipse.collections.impl.multimap.list.FastListMultimap;
import org.jetbrains.annotations.NotNull;

public class TestCiderGeneDatastore implements CiderGeneDatastore
{
    private final MutableMultimap<String, VJAnchorTemplate> mAnchorSequenceMap = new FastListMultimap<>();

    private final MutableMap<VJGeneType, MutableMultimap<String, VJAnchorTemplate>> mGeneTypeAnchorSeqMap = new UnifiedMap<>();
    private final MutableMultimap<VJAnchorReferenceLocation, VJAnchorTemplate> mGeneLocationVJGeneMap = new FastListMultimap<>();

    private final MutableList<IgTcrConstantRegion> mIgTcrConstantRegions = new FastList<>();

    @Override
    public SetIterable<String> getAnchorSequenceSet(@NotNull VJGeneType geneType)
    {
        Multimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.keySet() : Sets.immutable.empty();
    }

    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorSequence(@NotNull String anchorSeq)
    {
        return mAnchorSequenceMap.get(anchorSeq).toImmutableList();
    }

    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq)
    {
        Multimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.get(anchorSeq).toImmutableList() : Lists.immutable.empty();
    }

    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation)
    {
        return mGeneLocationVJGeneMap.get(vjAnchorReferenceLocation).toImmutableList();
    }

    @Override
    public SetIterable<VJAnchorReferenceLocation> getVJAnchorReferenceLocations()
    {
        return mGeneLocationVJGeneMap.keySet();
    }

    @Override
    public ImmutableCollection<IgTcrConstantRegion> getIgConstantRegions() { return mIgTcrConstantRegions.toImmutable(); }

    public TestCiderGeneDatastore(List<VJAnchorTemplate> vjAnchorTemplates)
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
                mGeneTypeAnchorSeqMap.computeIfAbsent(gene.getType(), o -> new FastListMultimap<>())
                        .put(gene.getAnchorSequence(), gene);
            }
        }
    }
}