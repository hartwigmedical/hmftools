package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class StructuralVariantImpliedCopyNumber {

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantImpliedCopyNumber.class);

    private final StructuralVariantPloidyFactory<PurpleCopyNumber> structuralVariantPloidyFactory;

    public StructuralVariantImpliedCopyNumber(final PurityAdjuster purityAdjuster) {
        this.structuralVariantPloidyFactory = new StructuralVariantPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
    }

    @NotNull
    public ListMultimap<String, PurpleCopyNumber> svImpliedCopyNumber(final List<StructuralVariant> structuralVariants,
            @NotNull final ListMultimap<String, PurpleCopyNumber> copyNumbers) {

        final ListMultimap<String, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(copyNumbers);

        long previousMissingCopyNumbers = result.size();
        long currentMissingCopyNumbers = missingCopyNumberCount(result);

        while (currentMissingCopyNumbers < previousMissingCopyNumbers && currentMissingCopyNumbers > 0) {
            final GenomePositionSelector<StructuralVariantPloidy> selector =
                    GenomePositionSelectorFactory.create(createPloidies(structuralVariants, result));

            for (Chromosome chromosome : HumanChromosome.values()) {
                final String chromosomeName = chromosome.toString();
                final List<PurpleCopyNumber> chromosomeCopyNumbers = result.get(chromosomeName);
                for (int i = 0; i < chromosomeCopyNumbers.size(); i++) {
                    final PurpleCopyNumber copyNumber = chromosomeCopyNumbers.get(i);
                    if (implyCopyNumberFromSV(copyNumber)) {
                        final Optional<StructuralVariantPloidy> optionalStart =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.start()));

                        final Optional<StructuralVariantPloidy> optionalEnd =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.end() + 1));

                        if (optionalStart.isPresent() || optionalEnd.isPresent()) {
                            chromosomeCopyNumbers.set(i, impliedCopyNumber(copyNumber, optionalStart, optionalEnd));
                        }
                    }
                }
            }

            previousMissingCopyNumbers = currentMissingCopyNumbers;
            currentMissingCopyNumbers = missingCopyNumberCount(result);
        }

        return result;
    }

    @VisibleForTesting
    @NotNull
    static PurpleCopyNumber impliedCopyNumber(final PurpleCopyNumber copyNumber, final Optional<StructuralVariantPloidy> start,
            final Optional<StructuralVariantPloidy> end) {
        final double startWeight = start.map(StructuralVariantPloidy::impliedRightCopyNumberWeight).orElse(0d);
        final double startCopyNumber = start.map(StructuralVariantPloidy::impliedRightCopyNumber).orElse(0d);

        final double endWeight = end.map(StructuralVariantPloidy::impliedLeftCopyNumberWeight).orElse(0d);
        final double endCopyNumber = end.map(StructuralVariantPloidy::impliedLeftCopyNumber).orElse(0d);

        final double newCopyNumber = (startCopyNumber * startWeight + endCopyNumber * endWeight) / (startWeight + endWeight);
        return ImmutablePurpleCopyNumber.builder().from(copyNumber).inferred(true).averageTumorCopyNumber(newCopyNumber).build();
    }

    @NotNull
    private List<StructuralVariantPloidy> createPloidies(final List<StructuralVariant> structuralVariants,
            @NotNull ListMultimap<String, PurpleCopyNumber> copyNumbers) {
        return structuralVariantPloidyFactory.create(structuralVariants, copyNumbers);
    }

    private long missingCopyNumberCount(Multimap<String, PurpleCopyNumber> copyNumbers) {
        return copyNumbers.values().stream().filter(this::implyCopyNumberFromSV).count();
    }

    private boolean implyCopyNumberFromSV(@NotNull final PurpleCopyNumber copyNumber) {
        return Doubles.isZero(copyNumber.averageTumorCopyNumber()) && copyNumber.inferred();
    }
}
