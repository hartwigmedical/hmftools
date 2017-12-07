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
    ListMultimap<String, CombinedRegion> svImpliedCopyNumber(final List<StructuralVariant> structuralVariants,
            @NotNull final ListMultimap<String, CombinedRegion> copyNumbers) {

        long previousMissingCopyNumbers = copyNumbers.size();
        long currentMissingCopyNumbers = missingCopyNumberCount(copyNumbers);

        while (currentMissingCopyNumbers < previousMissingCopyNumbers && currentMissingCopyNumbers > 0) {
            final GenomePositionSelector<StructuralVariantLegPloidy> selector =
                    GenomePositionSelectorFactory.create(createPloidies(structuralVariants, copyNumbers));

            for (Chromosome chromosome : HumanChromosome.values()) {
                final String chromosomeName = chromosome.toString();
                final List<CombinedRegion> chromosomeCopyNumbers = copyNumbers.get(chromosomeName);
                boolean svInferred = false;
                for (final CombinedRegion copyNumber : chromosomeCopyNumbers) {
                    if (implyCopyNumberFromSV(copyNumber)) {
                        final Optional<StructuralVariantLegPloidy> optionalStart =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.start()));

                        final Optional<StructuralVariantLegPloidy> optionalEnd =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.end() + 1));

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


    private static void inferCopyNumberFromStructuralVariants(@NotNull final CombinedRegion region, final Optional<StructuralVariantLegPloidy> start,
            final Optional<StructuralVariantLegPloidy> end) {
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
            @NotNull ListMultimap<String, CombinedRegion> copyNumbers) {
        return structuralVariantPloidyFactory.create(structuralVariants, copyNumbers);
    }

    private long missingCopyNumberCount(Multimap<String, CombinedRegion> copyNumbers) {
        return copyNumbers.values().stream().filter(this::implyCopyNumberFromSV).count();
    }

    private boolean implyCopyNumberFromSV(@NotNull final CombinedRegion copyNumber) {
        return !copyNumber.isProcessed();
    }
}
