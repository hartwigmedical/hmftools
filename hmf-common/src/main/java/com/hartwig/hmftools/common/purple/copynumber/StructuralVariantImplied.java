package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
class StructuralVariantImplied {

    private final StructuralVariantLegPloidyFactory<CombinedRegion> structuralVariantPloidyFactory;

    StructuralVariantImplied(final PurityAdjuster purityAdjuster) {
        this.structuralVariantPloidyFactory =
                new StructuralVariantLegPloidyFactory<>(purityAdjuster, x -> x.isProcessed() ? x.tumorCopyNumber() : 0);
    }

    @NotNull
    ListMultimap<Chromosome, CombinedRegion> svImpliedCopyNumber(final List<StructuralVariant> structuralVariants,
            @NotNull final ListMultimap<Chromosome, CombinedRegion> copyNumbers) {

        long previousMissingCopyNumbers = copyNumbers.size();
        long currentMissingCopyNumbers = missingCopyNumberCount(copyNumbers);

        while (currentMissingCopyNumbers < previousMissingCopyNumbers && currentMissingCopyNumbers > 0) {
            final GenomePositionSelector<StructuralVariantLegPloidy> selector =
                    GenomePositionSelectorFactory.create(createPloidies(structuralVariants, copyNumbers));

            for (Chromosome chromosome : HumanChromosome.values()) {
                final List<CombinedRegion> chromosomeCopyNumbers = copyNumbers.get(chromosome);
                boolean svInferred = false;
                for (final CombinedRegion copyNumber : chromosomeCopyNumbers) {
                    if (implyCopyNumberFromSV(copyNumber)) {
                        final Optional<StructuralVariantLegPloidy> optionalStart =
                                select(copyNumber.chromosome(), copyNumber.start(), selector);
                        final Optional<StructuralVariantLegPloidy> optionalEnd =
                                select(copyNumber.chromosome(), copyNumber.end() + 1, selector);
                        if (optionalStart.isPresent() || optionalEnd.isPresent()) {
                            svInferred = true;
                            inferCopyNumberFromStructuralVariants(copyNumber, optionalStart, optionalEnd);
                        }
                    }
                }

                // JOBA: Extend structural variant segments
                if (svInferred) {
                    ExtendStructuralVariant.extendStructuralVariants(chromosomeCopyNumbers);
                }
            }

            previousMissingCopyNumbers = currentMissingCopyNumbers;
            currentMissingCopyNumbers = missingCopyNumberCount(copyNumbers);
        }

        return copyNumbers;
    }

    @NotNull
    private Optional<StructuralVariantLegPloidy> select(@NotNull final String chromosome, long position,
            @NotNull final GenomePositionSelector<StructuralVariantLegPloidy> selector) {

        final Optional<StructuralVariantLegPloidy> posOrientation =
                selector.select(GenomePositions.create(chromosome, position - 1)).filter(x -> x.orientation() == 1);
        if (posOrientation.isPresent()) {
            return posOrientation;
        }

        return selector.select(GenomePositions.create(chromosome, position)).filter(x -> x.orientation() == -1);
    }

    private static void inferCopyNumberFromStructuralVariants(@NotNull final CombinedRegion region,
            final Optional<StructuralVariantLegPloidy> start, final Optional<StructuralVariantLegPloidy> end) {
        region.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, inferCopyNumberFromStructuralVariants(start, end));
    }

    @VisibleForTesting
    static double inferCopyNumberFromStructuralVariants(final Optional<StructuralVariantLegPloidy> start,
            final Optional<StructuralVariantLegPloidy> end) {
        final double startWeight = start.map(StructuralVariantLegPloidy::impliedRightCopyNumberWeight).orElse(0d);
        final double startCopyNumber = start.map(StructuralVariantLegPloidy::impliedRightCopyNumber).orElse(0d);

        final double endWeight = end.map(StructuralVariantLegPloidy::impliedLeftCopyNumberWeight).orElse(0d);
        final double endCopyNumber = end.map(StructuralVariantLegPloidy::impliedLeftCopyNumber).orElse(0d);

        return (startCopyNumber * startWeight + endCopyNumber * endWeight) / (startWeight + endWeight);
    }

    @NotNull
    private List<StructuralVariantLegPloidy> createPloidies(final List<StructuralVariant> structuralVariants,
            @NotNull ListMultimap<Chromosome, CombinedRegion> copyNumbers) {
        return structuralVariantPloidyFactory.create(structuralVariants, copyNumbers);
    }

    private long missingCopyNumberCount(Multimap<?, CombinedRegion> copyNumbers) {
        return copyNumbers.values().stream().filter(this::implyCopyNumberFromSV).count();
    }

    private boolean implyCopyNumberFromSV(@NotNull final CombinedRegion copyNumber) {
        return !copyNumber.isProcessed();
    }
}
