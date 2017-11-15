package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class StructuralVariantPloidyFactory {

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantPloidyFactory.class);

    @NotNull
    static List<StructuralVariantPloidy> create(@NotNull final StructuralVariant variant,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers) {

        if (variant.type() == StructuralVariantType.INS) {
            return Collections.emptyList();
        }

        final Optional<ModifiableStructuralVariantPloidy> start = Optional.ofNullable(variant.startAF())
                .flatMap(vaf -> create(variant.startChromosome(),
                        variant.startPosition(),
                        (int) variant.startOrientation(),
                        vaf,
                        GenomeRegionSelectorFactory.create(copyNumbers)));

        final Optional<ModifiableStructuralVariantPloidy> end = Optional.ofNullable(variant.endAF())
                .flatMap(vaf -> create(variant.endChromosome(),
                        variant.endPosition(),
                        (int) variant.endOrientation(),
                        vaf,
                        GenomeRegionSelectorFactory.create(copyNumbers)));

        if (!start.isPresent() && !end.isPresent()) {
            return Collections.emptyList();
        }

        final List<StructuralVariantPloidy> result = Lists.newArrayList();
        double startWeight = start.map(ModifiableStructuralVariantPloidy::weight).orElse(0D);
        double startPloidy = start.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D);
        double endWeight = end.map(ModifiableStructuralVariantPloidy::weight).orElse(0D);
        double endPloidy = end.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D);

        double totalWeight = startWeight + endWeight;
        double averagePloidy = (startWeight * startPloidy + endWeight * endPloidy) / totalWeight;

        if (start.isPresent()) {
            result.add(start.get().setWeight(totalWeight).setAverageImpliedPloidy(averagePloidy));
        } else {
            LOGGER.debug("Unable to determine start ploidy of {}:{}", variant.startChromosome(), variant.startPosition());
        }

        if (end.isPresent()) {
            result.add(end.get().setWeight(totalWeight).setAverageImpliedPloidy(averagePloidy));
        } else {
            LOGGER.debug("Unable to determine end ploidy of {}:{}", variant.startChromosome(), variant.startPosition());
        }

        Collections.sort(result);
        return result;
    }

    private static Optional<ModifiableStructuralVariantPloidy> create(@NotNull final String chromosome, long position, int orientation,
            double vaf, @NotNull final GenomeRegionSelector<PurpleCopyNumber> selector) {

        final GenomePosition svPosition = GenomePositions.create(chromosome, position);
        final GenomePosition svPositionLeft = GenomePositions.create(chromosome, position - 1);
        final Optional<PurpleCopyNumber> left = selector.select(svPositionLeft).filter(x -> !Doubles.isZero(x.averageTumorCopyNumber()));
        final Optional<PurpleCopyNumber> right = selector.select(svPosition).filter(x -> !Doubles.isZero(x.averageTumorCopyNumber()));

        final Optional<PurpleCopyNumber> correct;
        final Optional<PurpleCopyNumber> alternate;
        if (orientation == 1) {
            correct = left;
            alternate = right;
        } else {
            correct = right;
            alternate = left;
        }

        if (!correct.isPresent() && !alternate.isPresent()) {
            return Optional.empty();
        }

        final double ploidy;
        final double weight;
        if (correct.isPresent()) {
            ploidy = vaf * correct.get().averageTumorCopyNumber();
            weight = 1;
        } else {
            double copyNumber = alternate.get().averageTumorCopyNumber();
            ploidy = vaf / (1 - vaf) * copyNumber;
            weight = 1 / (1 + Math.pow(Math.max(copyNumber, 2) / Math.min(Math.max(copyNumber, 0.01), 2), 2));
        }

        return Optional.of(ModifiableStructuralVariantPloidy.create()
                .from(svPosition)
                .setOrientation(orientation)
                .setUnweightedImpliedPloidy(ploidy)
                .setLeftCopyNumber(left.map(PurpleCopyNumber::averageTumorCopyNumber))
                .setRightCopyNumber(right.map(PurpleCopyNumber::averageTumorCopyNumber))
                .setWeight(weight));
    }

}
