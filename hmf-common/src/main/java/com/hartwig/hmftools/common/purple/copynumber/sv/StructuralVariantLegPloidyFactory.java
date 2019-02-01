package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantLegPloidyFactory<T extends GenomeRegion> {

    private static final double VAF_TO_USE_READ_DEPTH = 0.75;
    private static final double RECIPROCAL_VAF_TO_USE_READ_DEPTH = reciprocalVAF(VAF_TO_USE_READ_DEPTH);

    private final double averageCopyNumber;
    private final int averageReadDepth;

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final StructuralVariantLegCopyNumberFactory<T> copyNumberFactory;

    public StructuralVariantLegPloidyFactory(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Function<T, Double> copyNumberExtractor) {
        this(0, 0, purityAdjuster, copyNumberExtractor);
    }

    public StructuralVariantLegPloidyFactory(int averageReadDepth, double averageCopyNumber, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Function<T, Double> copyNumberExtractor) {
        this.averageCopyNumber = averageCopyNumber;
        this.averageReadDepth = averageReadDepth;
        this.purityAdjuster = purityAdjuster;
        this.copyNumberFactory = new StructuralVariantLegCopyNumberFactory<>(copyNumberExtractor);
    }

    @NotNull
    public List<StructuralVariantLegPloidy> create(@NotNull final StructuralVariant variant,
            @NotNull final Multimap<Chromosome, T> copyNumbers) {
        final List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        final List<StructuralVariantLegs> allLegs = StructuralVariantLegsFactory.create(variant);

        for (StructuralVariantLegs leg : allLegs) {
            result.addAll(create(leg, copyNumbers));
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    public List<StructuralVariantLegPloidy> create(@NotNull final List<StructuralVariant> variants,
            @NotNull final Multimap<Chromosome, T> copyNumbers) {
        final List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        final List<StructuralVariantLegs> allLegs = StructuralVariantLegsFactory.create(variants);

        for (StructuralVariantLegs leg : allLegs) {
            result.addAll(create(leg, copyNumbers));
        }

        Collections.sort(result);
        return result;
    }

    @VisibleForTesting
    @NotNull
    List<StructuralVariantLegPloidy> create(@NotNull final StructuralVariantLegs legs, @NotNull final Multimap<Chromosome, T> copyNumbers) {
        final Optional<ModifiableStructuralVariantLegPloidy> start =
                legs.start().flatMap(x -> create(x, GenomeRegionSelectorFactory.createImproved(copyNumbers)));

        final Optional<ModifiableStructuralVariantLegPloidy> end =
                legs.end().flatMap(x -> create(x, GenomeRegionSelectorFactory.createImproved(copyNumbers)));

        if (!start.isPresent() && !end.isPresent()) {
            return Collections.emptyList();
        }

        final List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        double startWeight = start.map(ModifiableStructuralVariantLegPloidy::weight).orElse(0D);
        double startPloidy = start.map(ModifiableStructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D);
        double endWeight = end.map(ModifiableStructuralVariantLegPloidy::weight).orElse(0D);
        double endPloidy = end.map(ModifiableStructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D);

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
    Optional<ModifiableStructuralVariantLegPloidy> create(@NotNull final StructuralVariantLeg leg,
            @NotNull final GenomeRegionSelector<T> selector) {

        final StructuralVariantLegCopyNumber legCopyNumber = copyNumberFactory.create(leg, selector);

        final Optional<Double> largerCopyNumber;
        final Optional<Double> smallerCopyNumber;
        if (leg.orientation() == 1) {
            largerCopyNumber = legCopyNumber.leftCopyNumber();
            smallerCopyNumber = legCopyNumber.rightCopyNumber();
        } else {
            largerCopyNumber = legCopyNumber.rightCopyNumber();
            smallerCopyNumber = legCopyNumber.leftCopyNumber();
        }

        if (!largerCopyNumber.isPresent() && !smallerCopyNumber.isPresent()) {
            return Optional.empty();
        }

        final Double observedVaf = leg.alleleFrequency();
        if (observedVaf == null) {
            return Optional.empty();
        }

        final double adjustedVaf;
        final double ploidy;
        final double weight;

        if (largerCopyNumber.isPresent()) {
            double copyNumber = largerCopyNumber.get();
            adjustedVaf = purityAdjuster.purityAdjustedVAF(leg.chromosome(), Math.max(0.001, copyNumber), observedVaf);
            ploidy = adjustedVaf * copyNumber;
            weight = 1;
        } else {
            double copyNumber = smallerCopyNumber.get();
            double reciprocalVAF = reciprocalVAF(observedVaf);
            if (!Double.isFinite(reciprocalVAF)) {
                return Optional.empty();
            }

            adjustedVaf = purityAdjuster.purityAdjustedVAF(leg.chromosome(), Math.max(0.001, copyNumber), reciprocalVAF);

            if (averageReadDepth > 0 && Doubles.greaterThan(adjustedVaf, RECIPROCAL_VAF_TO_USE_READ_DEPTH)) {
                final Integer tumourVariantFragmentCount = leg.tumourVariantFragmentCount();
                if (tumourVariantFragmentCount != null && tumourVariantFragmentCount > 0) {
                    ploidy = readDepthImpliedPloidy(tumourVariantFragmentCount);
                } else {
                    return Optional.empty();
                }

            } else {
                ploidy = adjustedVaf * copyNumber;
            }

            weight = 1 / (1 + Math.pow(Math.max(copyNumber, 2) / Math.min(Math.max(copyNumber, 0.01), 2), 2));
        }

        return Optional.of(ModifiableStructuralVariantLegPloidy.create()
                .from(legCopyNumber)
                .from(leg)
                .setObservedVaf(observedVaf)
                .setAdjustedVaf(adjustedVaf)
                .setOrientation(leg.orientation())
                .setUnweightedImpliedPloidy(ploidy)
                .setWeight(weight));
    }

    private static double reciprocalVAF(double vaf) {
        return vaf / (1 - vaf);
    }

    private double readDepthImpliedPloidy(int tumourVariantFragmentCount) {
        return averageCopyNumber * tumourVariantFragmentCount / averageReadDepth;
    }
}
