package com.hartwig.hmftools.cdr3;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public interface VJGeneStore
{
    List<VJGene> getVJGenes();

    Set<String> getAnchorSequenceSet();
    Set<String> getAnchorSequenceSet(@NotNull VJGeneType geneType);
    Collection<VJGene> getByAnchorSequence(@NotNull String anchorSeq);
    Collection<VJGene> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq);
    Collection<VJGene> getByAnchorAminoAcidSequence(@NotNull VJGeneType geneType, @NotNull String anchorAminoAcidSeq);
    Collection<VJGene> getByAnchorGeneLocation(@NotNull GeneLocation geneLocation);
    Collection<GeneLocation> getVJGenomeRegions();
}
