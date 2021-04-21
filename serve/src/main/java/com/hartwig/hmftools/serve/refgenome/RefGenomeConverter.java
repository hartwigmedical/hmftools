package com.hartwig.hmftools.serve.refgenome;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

class RefGenomeConverter {

    @NotNull
    private final LiftOver liftOver;
    @NotNull
    private final Map<String, String> geneMap;
    @NotNull
    private final Map<String, String> chromosomeMap;

    public RefGenomeConverter(@NotNull final LiftOver liftOver, @NotNull final Map<String, String> geneMap,
            @NotNull final Map<String, String> chromosomeMap) {
        this.liftOver = liftOver;
        this.geneMap = geneMap;
        this.chromosomeMap = chromosomeMap;
    }

    @NotNull
    public Set<KnownHotspot> convertKnownHotspots(@NotNull Set<KnownHotspot> hotspots) {
        Set<KnownHotspot> convertedHotspots = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            Interval interval = new Interval(RefGenomeFunctions.enforceChromosome(hotspot.chromosome()),
                    (int) hotspot.position(),
                    (int) hotspot.position());
            Interval lifted = liftOver.liftOver(interval);

            convertedHotspots.add(ImmutableKnownHotspot.builder()
                    .from(hotspot)
                    .chromosome(lifted.getContig())
                    .position(lifted.getStart())
                    .build());
        }
        return convertedHotspots;
    }

    @NotNull
    public Set<KnownCodon> convertKnownCodons(@NotNull Set<KnownCodon> codons) {
        return codons;
    }

    @NotNull
    public Set<KnownExon> convertKnownExons(@NotNull Set<KnownExon> exons) {
        return exons;
    }

    @NotNull
    public Set<KnownCopyNumber> convertKnownCopyNumbers(@NotNull Set<KnownCopyNumber> copyNumbers) {
        return copyNumbers;
    }

    @NotNull
    public Set<KnownFusionPair> convertKnownFusionPairs(@NotNull Set<KnownFusionPair> knownFusionPairs) {
        return knownFusionPairs;
    }

    @NotNull
    public Set<ActionableHotspot> convertActionableHotspots(@NotNull Set<ActionableHotspot> actionableHotspots) {
        return actionableHotspots;
    }

    @NotNull
    public Set<ActionableRange> convertActionableRanges(@NotNull Set<ActionableRange> actionableRanges) {
        return actionableRanges;
    }

    @NotNull
    public Set<ActionableGene> convertActionableGenes(@NotNull Set<ActionableGene> actionableGenes) {
        return actionableGenes;
    }

    @NotNull
    public Set<ActionableFusion> convertActionableFusion(@NotNull Set<ActionableFusion> actionableFusions) {
        return actionableFusions;
    }
}
