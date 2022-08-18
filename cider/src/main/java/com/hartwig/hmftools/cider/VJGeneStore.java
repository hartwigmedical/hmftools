package com.hartwig.hmftools.cider;

import java.util.Collection;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

// maybe rename to IG TCR template store
public interface VJGeneStore
{
    Set<String> getAnchorSequenceSet();
    Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType);
    Collection<VJAnchorTemplate> getByAnchorSequence(@NotNull String anchorSeq);
    Collection<VJAnchorTemplate> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq);
    Collection<VJAnchorTemplate> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation);
    Collection<VJAnchorReferenceLocation> getVJAnchorReferenceLocations();

    Collection<IgTcrConstantRegion> getIgConstantRegions();
}
