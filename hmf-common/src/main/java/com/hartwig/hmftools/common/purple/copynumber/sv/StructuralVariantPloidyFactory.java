package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantPloidyFactory<T extends GenomeRegion> {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final Function<T, Double> copyNumberExtractor;

    public StructuralVariantPloidyFactory(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Function<T, Double> copyNumberExtractor) {
        this.purityAdjuster = purityAdjuster;
        this.copyNumberExtractor = copyNumberExtractor;
    }

    @NotNull
    public List<StructuralVariantPloidy> create(@NotNull final List<StructuralVariant> variants,
            @NotNull final Multimap<String, T> copyNumbers) {
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
    List<StructuralVariantPloidy> create(@NotNull final StructuralVariantLegs legs, @NotNull final Multimap<String, T> copyNumbers) {
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

    @VisibleForTesting
    @NotNull
    Optional<ModifiableStructuralVariantPloidy> create(@NotNull StructuralVariantLeg leg, @NotNull final GenomeRegionSelector<T> selector) {
        final GenomePosition svPositionLeft = GenomePositions.create(leg.chromosome(), leg.position() - 1);
        final Optional<Double> left =
                selector.select(svPositionLeft).flatMap(x -> Optional.ofNullable(copyNumberExtractor.apply(x))).filter(Doubles::positive);
        final Optional<Double> right =
                selector.select(leg).flatMap(x -> Optional.ofNullable(copyNumberExtractor.apply(x))).filter(Doubles::positive);

        final Optional<Double> correct;
        final Optional<Double> alternate;
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

        final double vaf = leg.alleleFrequency();
        final double ploidy;
        final double weight;
        if (correct.isPresent()) {
            double copyNumber = correct.get();
            ploidy = purityAdjustedPloidy(leg.chromosome(), vaf, copyNumber);
            weight = 1;
        } else {
            double copyNumber = alternate.get();
            ploidy = purityAdjustedPloidy(leg.chromosome(), vaf / (1 - vaf), copyNumber);
            weight = 1 / (1 + Math.pow(Math.max(copyNumber, 2) / Math.min(Math.max(copyNumber, 0.01), 2), 2));
        }

        return Optional.of(ModifiableStructuralVariantPloidy.create()
                .from(leg)
                .setVaf(vaf)
                .setOrientation(leg.orientation())
                .setUnweightedImpliedPloidy(ploidy)
                .setLeftCopyNumber(left)
                .setRightCopyNumber(right)
                .setWeight(weight));
    }

    private double purityAdjustedPloidy(@NotNull String chromosome, double vaf, double copyNumber) {
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(chromosome, copyNumber, vaf);
        return adjustedVAF * copyNumber;
    }
}
