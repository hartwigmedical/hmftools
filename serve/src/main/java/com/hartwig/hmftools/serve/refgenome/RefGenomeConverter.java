package com.hartwig.hmftools.serve.refgenome;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.common.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.common.serve.datamodel.MutationTypeFilter;
import com.hartwig.hmftools.common.serve.datamodel.range.RangeAnnotation;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
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

    public RefGenomeConverter(@NotNull final RefGenomeVersion sourceVersion, @NotNull final RefGenomeVersion targetVersion,
            @NotNull final IndexedFastaSequenceFile targetSequence, @NotNull final LiftOverAlgo liftOverAlgo) {
        this.sourceVersion = sourceVersion;
        this.targetVersion = targetVersion;
        this.targetSequence = targetSequence;
        this.liftOverAlgo = liftOverAlgo;
    }

    @NotNull
    public Set<KnownHotspot> convertKnownHotspots(@NotNull Set<KnownHotspot> hotspots) {
        Set<KnownHotspot> convertedHotspots = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            VariantHotspot lifted = liftOverHotspot(hotspot);

            if (lifted != null) {
                convertedHotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
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
                    // We blank out the transcript and codon rank since we are unsure to what extend
                    // the transcript maps to the new ref genome.
                    convertedCodons.add(ImmutableKnownCodon.builder()
                            .from(codon)
                            .annotation(ImmutableCodonAnnotation.builder().from(liftedAnnotation).build())
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
                // We blank out the transcript and exon rank since we are unsure to what extend the transcript maps to the new ref genome.
                convertedExons.add(ImmutableKnownExon.builder()
                        .from(exon)
                        .annotation(ImmutableExonAnnotation.builder().from(liftedAnnotation).build())
                        .build());
            }
        }

        return convertedExons;
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
                        .chromosome(lifted.chromosome())
                        .start(lifted.start())
                        .end(lifted.end())
                        .build());
            }
        }
        return convertedActionableRanges;
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

        // We blank out the transcript and rank since we are unsure to what extend the transcript maps to the new ref genome.
        return new RangeAnnotation() {
            @NotNull
            @Override
            public String gene() {
                return annotation.gene();
            }

            @NotNull
            @Override
            public String transcript() {
                return Strings.EMPTY;
            }

            @Override
            public int rank() {
                return 0;
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
            public int start() {
                return liftedStart.position();
            }

            @Override
            public int end() {
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
    private String sequence(@NotNull String chromosome, long start, long length) {
        String targetChromosome = targetVersion.versionedChromosome(chromosome);
        return targetSequence.getSubsequenceAt(targetChromosome, start, start + length - 1).getBaseString();
    }
}