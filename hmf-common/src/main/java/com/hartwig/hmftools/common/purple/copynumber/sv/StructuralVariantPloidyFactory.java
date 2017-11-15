package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

class StructuralVariantPloidyFactory {

    @NotNull
    static List<StructuralVariantPloidy> create(@NotNull final List<StructuralVariant> variants,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers) {

        final List<StructuralVariantPloidy> result = Lists.newArrayList();
        final List<StructuralVariantLegs> allLegs = StructuralVariantLegsFactory.create(variants);

        for (StructuralVariantLegs leg : allLegs) {
            result.addAll(create(leg, copyNumbers));
        }

        Collections.sort(result);
        return result;
    }

    @VisibleForTesting
    @NotNull
    static List<StructuralVariantPloidy> create(@NotNull final StructuralVariantLegs legs,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers) {

        final Optional<ModifiableStructuralVariantPloidy> start =
                legs.start().flatMap(x -> create(x, GenomeRegionSelectorFactory.create(copyNumbers)));

        final Optional<ModifiableStructuralVariantPloidy> end =
                legs.end().flatMap(x -> create(x, GenomeRegionSelectorFactory.create(copyNumbers)));

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

        start.ifPresent(modifiableStructuralVariantPloidy -> result.add(modifiableStructuralVariantPloidy.setWeight(totalWeight)
                .setAverageImpliedPloidy(averagePloidy)));

        end.ifPresent(modifiableStructuralVariantPloidy -> result.add(modifiableStructuralVariantPloidy.setWeight(totalWeight)
                .setAverageImpliedPloidy(averagePloidy)));

        Collections.sort(result);
        return result;
    }

    private static Optional<ModifiableStructuralVariantPloidy> create(@NotNull StructuralVariantLeg leg,
            @NotNull final GenomeRegionSelector<PurpleCopyNumber> selector) {

        final GenomePosition svPositionLeft = GenomePositions.create(leg.chromosome(), leg.position() - 1);
        final Optional<PurpleCopyNumber> left = selector.select(svPositionLeft).filter(x -> !Doubles.isZero(x.averageTumorCopyNumber()));
        final Optional<PurpleCopyNumber> right = selector.select(leg).filter(x -> !Doubles.isZero(x.averageTumorCopyNumber()));

        final Optional<PurpleCopyNumber> correct;
        final Optional<PurpleCopyNumber> alternate;
        if (leg.orientation() == 1) {
            correct = left;
            alternate = right;
        } else {
            correct = right;
            alternate = left;
        }

        if (!correct.isPresent() && !alternate.isPresent()) {
            return Optional.empty();
        }

        final double vaf = leg.vaf();
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
                .from(leg)
                .setOrientation(leg.orientation())
                .setUnweightedImpliedPloidy(ploidy)
                .setLeftCopyNumber(left.map(PurpleCopyNumber::averageTumorCopyNumber))
                .setRightCopyNumber(right.map(PurpleCopyNumber::averageTumorCopyNumber))
                .setWeight(weight));
    }
}
