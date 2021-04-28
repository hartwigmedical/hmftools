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
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;

class RefGenomeConverter {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeConverter.class);

    @NotNull
    private final RefGenomeVersion sourceVersion;
    @NotNull
    private final RefGenomeVersion targetVersion;
    @NotNull
    private final IndexedFastaSequenceFile targetSequence;
    @NotNull
    private final LiftOver liftOver;
    @NotNull
    private final GeneNameMapping geneNameMapping;

    public RefGenomeConverter(@NotNull final RefGenomeVersion sourceVersion, @NotNull final RefGenomeVersion targetVersion,
            @NotNull final IndexedFastaSequenceFile targetSequence, @NotNull final LiftOver liftOver,
            @NotNull final GeneNameMapping geneNameMapping) {
        this.sourceVersion = sourceVersion;
        this.targetVersion = targetVersion;
        this.targetSequence = targetSequence;
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

            if (lifted == null) {
                LOGGER.warn("Liftover could not be performed for '{}' on '{}'", hotspot.proteinAnnotation(), hotspot.gene());
            } else {
                String versionedOldChromosome = targetVersion.versionedChromosome(hotspot.chromosome());
                if (!lifted.getContig().equals(versionedOldChromosome)) {
                    LOGGER.warn("Liftover moved chromosome from '{}' to '{}' on {}", versionedOldChromosome, lifted.getContig(), hotspot);
                }

                String newRef = targetSequence.getSubsequenceAt(lifted.getContig(),
                        lifted.getStart(),
                        lifted.getStart() + hotspot.ref().length() - 1).getBaseString();
                if (!newRef.equals(hotspot.ref())) {
                    LOGGER.warn("Skipping liftover: Ref changed from '{}' to '{}' on {}", hotspot.ref(), newRef, hotspot);
                } else {
                    convertedHotspots.add(ImmutableKnownHotspot.builder()
                            .from(hotspot)
                            .gene(mapGene(hotspot.gene(), sourceVersion, targetVersion))
                            .chromosome(lifted.getContig())
                            .position(lifted.getStart())
                            .build());
                }
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
        String mappedGene;
        if (sourceVersion == targetVersion) {
            mappedGene = gene;
        } else if (sourceVersion == RefGenomeVersion.V37 && targetVersion == RefGenomeVersion.V38) {
            mappedGene = geneNameMapping.v38Gene(gene);
        } else if (sourceVersion == RefGenomeVersion.V38 && targetVersion == RefGenomeVersion.V37) {
            mappedGene = geneNameMapping.v37Gene(gene);
        } else {
            throw new IllegalStateException("Cannot map genes from ref genome version " + sourceVersion + " to " + targetVersion);
        }

        if (!mappedGene.equals(gene)) {
            LOGGER.info("Mapped gene '{}' for {} to '{}' on {}", gene, sourceVersion, mappedGene, targetVersion);
        }

        return mappedGene;
    }
}
