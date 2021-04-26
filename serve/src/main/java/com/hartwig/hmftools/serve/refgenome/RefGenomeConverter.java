package com.hartwig.hmftools.serve.refgenome;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

class RefGenomeConverter {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeConverter.class);

    @NotNull
    private final RefGenomeVersion sourceVersion;
    @NotNull
    private final RefGenomeVersion targetVersion;
    @NotNull
    private final LiftOver liftOver;
    @NotNull
    private final GeneNameMapping geneNameMapping;

    public RefGenomeConverter(@NotNull final RefGenomeVersion sourceVersion, @NotNull final RefGenomeVersion targetVersion,
            @NotNull final LiftOver liftOver, @NotNull final GeneNameMapping geneNameMapping) {
        assert sourceVersion != targetVersion;

        this.sourceVersion = sourceVersion;
        this.targetVersion = targetVersion;
        this.liftOver = liftOver;
        this.geneNameMapping = geneNameMapping;
    }

    @NotNull
    public Set<KnownHotspot> convertKnownHotspots(@NotNull Set<KnownHotspot> hotspots) {
        Set<KnownHotspot> convertedHotspots = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            Interval interval = new Interval(RefGenomeFunctions.enforceChromosome(hotspot.chromosome()),
                    (int) hotspot.position(),
                    (int) hotspot.position());
            Interval lifted = liftOver.liftOver(interval);

            if (lifted == null || lifted.getContig() == null) {
                LOGGER.warn("Liftover of '{}' led to non-interpretable '{}'", hotspot, lifted);
            }
            else {
                if (!sourceVersion.versionedChromosome(lifted.getContig()).equals(hotspot.chromosome())) {
                    LOGGER.warn("Liftover moved chromosome from '{}' to '{}' on {}", lifted.getContig(), hotspot.chromosome(), hotspot);
                }

                // TODO: Check if ref is the same.
                convertedHotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
                        .gene(mapGene(hotspot.gene(), sourceVersion, targetVersion))
                        .chromosome(targetVersion.versionedChromosome(hotspot.chromosome()))
                        .position(lifted.getStart())
                        .build());
            }
        }
        return convertedHotspots;
    }

    @NotNull
    public Set<KnownCodon> convertKnownCodons(@NotNull Set<KnownCodon> codons) {
        // TODO Implement
        return codons;
    }

    @NotNull
    public Set<KnownExon> convertKnownExons(@NotNull Set<KnownExon> exons) {
        // TODO Implement
        return exons;
    }

    @NotNull
    public Set<KnownCopyNumber> convertKnownCopyNumbers(@NotNull Set<KnownCopyNumber> copyNumbers) {
        // TODO Implement
        return copyNumbers;
    }

    @NotNull
    public Set<KnownFusionPair> convertKnownFusionPairs(@NotNull Set<KnownFusionPair> knownFusionPairs) {
        // TODO Implement
        return knownFusionPairs;
    }

    @NotNull
    public Set<ActionableHotspot> convertActionableHotspots(@NotNull Set<ActionableHotspot> actionableHotspots) {
        // TODO Implement
        return actionableHotspots;
    }

    @NotNull
    public Set<ActionableRange> convertActionableRanges(@NotNull Set<ActionableRange> actionableRanges) {
        // TODO Implement
        return actionableRanges;
    }

    @NotNull
    public Set<ActionableGene> convertActionableGenes(@NotNull Set<ActionableGene> actionableGenes) {
        // TODO Implement
        return actionableGenes;
    }

    @NotNull
    public Set<ActionableFusion> convertActionableFusion(@NotNull Set<ActionableFusion> actionableFusions) {
        // TODO Implement
        return actionableFusions;
    }

    @NotNull
    private String mapGene(@NotNull String gene, @NotNull RefGenomeVersion sourceVersion, @NotNull RefGenomeVersion targetVersion) {
        if (sourceVersion == targetVersion) {
            return gene;
        } else if (sourceVersion == RefGenomeVersion.V37 && targetVersion == RefGenomeVersion.V38) {
            return geneNameMapping.v38Gene(gene);
        } else if (sourceVersion == RefGenomeVersion.V38 && targetVersion == RefGenomeVersion.V37) {
            return geneNameMapping.v37Gene(gene);
        } else {
            throw new IllegalStateException("Cannot map genes from ref genome version " + sourceVersion + " to " + targetVersion);
        }
    }
}
