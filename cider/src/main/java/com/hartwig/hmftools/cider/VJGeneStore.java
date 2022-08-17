package com.hartwig.hmftools.cider;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

// maybe rename to IG locus data store?
public interface VJGeneStore
{
    List<VJGene> getVJGenes();

    Set<String> getAnchorSequenceSet();
    Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType);
    Collection<VJGene> getByAnchorSequence(@NotNull String anchorSeq);
    Collection<VJGene> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq);
    Collection<VJGene> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation);
    Collection<VJAnchorReferenceLocation> getVJAnchorReferenceLocations();

    Collection<IgConstantRegion> getIgConstantRegions();
}
