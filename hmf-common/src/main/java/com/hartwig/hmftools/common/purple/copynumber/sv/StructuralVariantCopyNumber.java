package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class StructuralVariantCopyNumber {

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantCopyNumber.class);

    private final List<StructuralVariant> structuralVariants;

    public StructuralVariantCopyNumber(final List<StructuralVariant> structuralVariants) {
        this.structuralVariants = structuralVariants;
    }

    @NotNull
    public ListMultimap<String, PurpleCopyNumber> calculateSVCopyNumber(@NotNull final ListMultimap<String, PurpleCopyNumber> copyNumbers) {

        final ListMultimap<String, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(copyNumbers);

        long previousMissingCopyNumbers = result.size();
        long currentMissingCopyNumbers = missingCopyNumbers(result);

        while (currentMissingCopyNumbers < previousMissingCopyNumbers && currentMissingCopyNumbers > 0) {
            final GenomePositionSelector<StructuralVariantPloidy> selector = GenomePositionSelectorFactory.create(createPloidies(result));

            for (Chromosome chromosome : HumanChromosome.values()) {
                final String chromosomeName = chromosome.toString();
                final List<PurpleCopyNumber> chromosomeCopyNumbers = result.get(chromosomeName);
                for (int i = 0; i < chromosomeCopyNumbers.size(); i++) {
                    final PurpleCopyNumber copyNumber = chromosomeCopyNumbers.get(i);
                    if (Doubles.isZero(copyNumber.averageTumorCopyNumber())) {
                        final Optional<StructuralVariantPloidy> optionalStart =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.start()));

                        final Optional<StructuralVariantPloidy> optionalEnd =
                                selector.select(GenomePositions.create(chromosomeName, copyNumber.end() + 1));

                        if (optionalStart.isPresent() || optionalEnd.isPresent()) {

                            final double startWeight = optionalStart.map(StructuralVariantPloidy::impliedRightCopyNumberWeight).orElse(0d);
                            final double startCopyNumber = optionalStart.map(StructuralVariantPloidy::impliedRightCopyNumber).orElse(0d);

                            final double endWeight = optionalEnd.map(StructuralVariantPloidy::impliedLeftCopyNumberWeight).orElse(0d);
                            final double endCopyNumber = optionalEnd.map(StructuralVariantPloidy::impliedLeftCopyNumber).orElse(0d);

                            final double newCopyNumber =
                                    (startCopyNumber * startWeight + endCopyNumber * endWeight) / (startWeight + endWeight);
                            final PurpleCopyNumber newPurpleCopyNumber = ImmutablePurpleCopyNumber.builder()
                                    .from(copyNumber)
                                    .inferred(true)
                                    .averageTumorCopyNumber(newCopyNumber)
                                    .build();

                            chromosomeCopyNumbers.set(i, newPurpleCopyNumber);
                        }

                    }

                }
            }

            previousMissingCopyNumbers = currentMissingCopyNumbers;
            currentMissingCopyNumbers = missingCopyNumbers(result);
        }

        return result;
    }

    @NotNull
    private List<StructuralVariantPloidy> createPloidies(@NotNull ListMultimap<String, PurpleCopyNumber> copyNumbers) {
        return StructuralVariantPloidyFactory.create(structuralVariants, copyNumbers);
    }

    private long missingCopyNumbers(Multimap<String, PurpleCopyNumber> copyNumbers) {
        return copyNumbers.values().stream().filter(x -> Doubles.isZero(x.averageTumorCopyNumber())).count();
    }

}
