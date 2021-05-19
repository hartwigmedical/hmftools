package com.hartwig.hmftools.serve.refgenome;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.extraction.range.RangeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverAlgo;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverChecker;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class RefGenomeConverter {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeConverter.class);

    @NotNull
    private final RefGenomeVersion sourceVersion;
    @NotNull
    private final RefGenomeVersion targetVersion;
    @NotNull
    private final IndexedFastaSequenceFile targetSequence;
    @NotNull
    private final LiftOverAlgo liftOverAlgo;
    @NotNull
    private final GeneNameMapping geneNameMapping;

    public RefGenomeConverter(@NotNull final RefGenomeVersion sourceVersion, @NotNull final RefGenomeVersion targetVersion,
            @NotNull final IndexedFastaSequenceFile targetSequence, @NotNull final LiftOverAlgo liftOverAlgo,
            @NotNull final GeneNameMapping geneNameMapping) {
        this.sourceVersion = sourceVersion;
        this.targetVersion = targetVersion;
        this.targetSequence = targetSequence;
        this.liftOverAlgo = liftOverAlgo;
        this.geneNameMapping = geneNameMapping;
    }

    @NotNull
    public Set<KnownHotspot> convertKnownHotspots(@NotNull Set<KnownHotspot> hotspots) {
        Set<KnownHotspot> convertedHotspots = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            VariantHotspot lifted = liftOverHotspot(hotspot);

            if (lifted != null) {
                convertedHotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
                        .gene(mapGene(hotspot.gene()))
                        .chromosome(lifted.chromosome())
                        .position(lifted.position())
                        .build());
            }
        }

        return convertedHotspots;
    }

    @NotNull
    public Set<KnownCodon> convertKnownCodons(@NotNull Set<KnownCodon> codons) {
        Set<KnownCodon> convertedCodons = Sets.newHashSet();
        for (KnownCodon codon : codons) {
            RangeAnnotation originalAnnotation = codon.annotation();
            RangeAnnotation liftedAnnotation = liftOverRange(originalAnnotation);
            if (liftedAnnotation != null) {
                if (originalAnnotation.end() - originalAnnotation.start() == 2 && liftedAnnotation.end() - liftedAnnotation.start() != 2) {
                    LOGGER.warn(" Skipping liftover from {} to {}: Lifted codon '{}' is no longer 3 bases long. Lifted codon: '{}'",
                            sourceVersion,
                            targetVersion,
                            originalAnnotation,
                            liftedAnnotation);
                } else {
                    // We blank out the transcript and codon index since we are unsure to what extend
                    // the transcript maps to the new ref genome.
                    convertedCodons.add(ImmutableKnownCodon.builder()
                            .from(codon)
                            .annotation(ImmutableCodonAnnotation.builder()
                                    .from(liftedAnnotation)
                                    .transcript(Strings.EMPTY)
                                    .codonIndex(0)
                                    .build())
                            .build());
                }
            }
        }

        return convertedCodons;
    }

    @NotNull
    public Set<KnownExon> convertKnownExons(@NotNull Set<KnownExon> exons) {
        Set<KnownExon> convertedExons = Sets.newHashSet();
        for (KnownExon exon : exons) {
            RangeAnnotation liftedAnnotation = liftOverRange(exon.annotation());
            if (liftedAnnotation != null) {
                // We blank out the transcript and exon index since we are unsure to what extend the transcript maps to the new ref genome.
                convertedExons.add(ImmutableKnownExon.builder()
                        .from(exon)
                        .annotation(ImmutableExonAnnotation.builder().from(liftedAnnotation).transcript(Strings.EMPTY).exonIndex(0).build())
                        .build());
            }
        }

        return convertedExons;
    }

    @NotNull
    public Set<KnownCopyNumber> convertKnownCopyNumbers(@NotNull Set<KnownCopyNumber> copyNumbers) {
        Set<KnownCopyNumber> convertedCopyNumbers = Sets.newHashSet();
        for (KnownCopyNumber copyNumber : copyNumbers) {
            convertedCopyNumbers.add(ImmutableKnownCopyNumber.builder().from(copyNumber).gene(mapGene(copyNumber.gene())).build());
        }

        return convertedCopyNumbers;
    }

    @NotNull
    public Set<KnownFusionPair> convertKnownFusionPairs(@NotNull Set<KnownFusionPair> fusionPairs) {
        Set<KnownFusionPair> convertedFusionPairs = Sets.newHashSet();
        for (KnownFusionPair fusionPair : fusionPairs) {
            convertedFusionPairs.add(ImmutableKnownFusionPair.builder()
                    .from(fusionPair)
                    .geneUp(mapGene(fusionPair.geneUp()))
                    .geneDown(mapGene(fusionPair.geneDown()))
                    .build());
        }

        return convertedFusionPairs;
    }

    @NotNull
    public Set<ActionableHotspot> convertActionableHotspots(@NotNull Set<ActionableHotspot> actionableHotspots) {
        Set<ActionableHotspot> convertedActionableHotspots = Sets.newHashSet();
        for (ActionableHotspot actionableHotspot : actionableHotspots) {
            VariantHotspot lifted = liftOverHotspot(actionableHotspot);
            if (lifted != null) {
                convertedActionableHotspots.add(ImmutableActionableHotspot.builder()
                        .from(actionableHotspot)
                        .chromosome(lifted.chromosome())
                        .position(lifted.position())
                        .build());
            }
        }
        return convertedActionableHotspots;
    }

    @NotNull
    public Set<ActionableRange> convertActionableRanges(@NotNull Set<ActionableRange> actionableRanges) {
        Set<ActionableRange> convertedActionableRanges = Sets.newHashSet();
        for (ActionableRange actionableRange : actionableRanges) {
            RangeAnnotation lifted = liftOverRange(actionableRange);
            if (lifted != null) {
                convertedActionableRanges.add(ImmutableActionableRange.builder()
                        .from(actionableRange)
                        .gene(lifted.gene())
                        .chromosome(lifted.chromosome())
                        .start(lifted.start())
                        .end(lifted.end())
                        .build());
            }
        }
        return convertedActionableRanges;
    }

    @NotNull
    public Set<ActionableGene> convertActionableGenes(@NotNull Set<ActionableGene> actionableGenes) {
        Set<ActionableGene> convertedActionableGenes = Sets.newHashSet();
        for (ActionableGene actionableGene : actionableGenes) {
            convertedActionableGenes.add(ImmutableActionableGene.builder()
                    .from(actionableGene)
                    .gene(mapGene(actionableGene.gene()))
                    .build());
        }
        return convertedActionableGenes;
    }

    @NotNull
    public Set<ActionableFusion> convertActionableFusions(@NotNull Set<ActionableFusion> actionableFusions) {
        Set<ActionableFusion> convertedActionableFusions = Sets.newHashSet();
        for (ActionableFusion actionableFusion : actionableFusions) {
            convertedActionableFusions.add(ImmutableActionableFusion.builder()
                    .from(actionableFusion)
                    .geneUp(mapGene(actionableFusion.geneUp()))
                    .geneDown(mapGene(actionableFusion.geneDown()))
                    .build());
        }
        return convertedActionableFusions;
    }

    @Nullable
    private VariantHotspot liftOverHotspot(@NotNull VariantHotspot hotspot) {
        LiftOverResult lifted = liftOverAlgo.liftOver(hotspot.chromosome(), hotspot.position());

        if (!LiftOverChecker.isValidLiftedPosition(lifted, hotspot)) {
            return null;
        }

        verifyNoChromosomeChange(hotspot.chromosome(), lifted, hotspot);

        String newRef = sequence(lifted.chromosome(), lifted.position(), hotspot.ref().length());
        if (!newRef.equals(hotspot.ref())) {
            LOGGER.warn(" Skipping liftover from {} to {}: Ref changed from '{}' to '{}' on position {} from {}",
                    sourceVersion,
                    targetVersion,
                    hotspot.ref(),
                    newRef,
                    lifted.position(),
                    hotspot);
            return null;
        }

        return ImmutableVariantHotspotImpl.builder().from(hotspot).chromosome(lifted.chromosome()).position(lifted.position()).build();
    }

    @Nullable
    private RangeAnnotation liftOverRange(@NotNull RangeAnnotation annotation) {
        LiftOverResult liftedStart = liftOverAlgo.liftOver(annotation.chromosome(), annotation.start());
        LiftOverResult liftedEnd = liftOverAlgo.liftOver(annotation.chromosome(), annotation.end());

        if (!LiftOverChecker.isValidLiftedRegion(liftedStart, liftedEnd, annotation)) {
            return null;
        }

        verifyNoChromosomeChange(annotation.chromosome(), liftedStart, annotation);
        verifyNoChromosomeChange(annotation.chromosome(), liftedEnd, annotation);

        // We blank out the transcript since we are unsure to what extend the transcript maps to the new ref genome.
        return new RangeAnnotation() {
            @NotNull
            @Override
            public String gene() {
                return mapGene(annotation.gene());
            }

            @NotNull
            @Override
            public MutationTypeFilter mutationType() {
                return annotation.mutationType();
            }

            @NotNull
            @Override
            public String chromosome() {
                return liftedStart.chromosome();
            }

            @Override
            public long start() {
                return liftedStart.position();
            }

            @Override
            public long end() {
                return liftedEnd.position();
            }

            @Override
            public String toString() {
                return chromosome() + ":" + start() + "-" + end();
            }
        };
    }

    private void verifyNoChromosomeChange(@NotNull String prevChromosome, @NotNull LiftOverResult lifted, @NotNull Object object) {
        String versionedChromosome = targetVersion.versionedChromosome(prevChromosome);
        if (!lifted.chromosome().equals(versionedChromosome)) {
            LOGGER.warn(" Liftover from {} to {} moved chromosome from '{}' to '{}' on {}",
                    sourceVersion,
                    targetVersion,
                    versionedChromosome,
                    lifted.chromosome(),
                    object);
        }
    }

    @NotNull
    private String mapGene(@NotNull String gene) {
        String mappedGene;
        if (sourceVersion == targetVersion) {
            mappedGene = gene;
        } else if (sourceVersion == RefGenomeVersion.V37 && targetVersion == RefGenomeVersion.V38) {
            if (geneNameMapping.isValidV37Gene(gene)) {
                mappedGene = geneNameMapping.v38Gene(gene);
            } else {
                LOGGER.warn(" Not a valid 37 gene encountered during liftover: '{}'!", gene);
                mappedGene = gene;
            }
        } else if (sourceVersion == RefGenomeVersion.V38 && targetVersion == RefGenomeVersion.V37) {
            if (geneNameMapping.isValidV38Gene(gene)) {
                mappedGene = geneNameMapping.v37Gene(gene);
            } else {
                LOGGER.warn(" Not a valid 38 gene encountered during liftover: '{}'!", gene);
                mappedGene = gene;
            }
        } else {
            throw new IllegalStateException("Cannot map genes from ref genome version " + sourceVersion + " to " + targetVersion);
        }

        if (!mappedGene.equals(gene)) {
            LOGGER.debug(" Mapped gene '{}' for {} to '{}' on {}", gene, sourceVersion, mappedGene, targetVersion);
            if (mappedGene.equals("NA")) {
                LOGGER.warn(" Gene '{}' mapped to 'NA' from {} to {}", gene, sourceVersion, targetVersion);
            }
        }

        return mappedGene;
    }

    @NotNull
    private String sequence(@NotNull String chromosome, long start, long length) {
        String targetChromosome = targetVersion.versionedChromosome(chromosome);
        return targetSequence.getSubsequenceAt(targetChromosome, start, start + length - 1).getBaseString();
    }
}
