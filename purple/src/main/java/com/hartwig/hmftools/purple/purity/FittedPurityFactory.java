package com.hartwig.hmftools.purple.purity;

import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterAbnormalChromosome;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Downsample;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.purple.region.FittedRegionFactory;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory {

    private static final int MAX_SOMATICS_TO_FIT = 1000;
    private static final double MAX_TUMOR_RATIO_TO_FIT = 3;

    private final double minPurity;
    private final double maxPurity;
    private final double purityIncrements;
    private final double somaticPenaltyWeight;
    private final CobaltChromosomes cobaltChromosomes;

    private final int totalBAFCount;
    private final double averageFittingRatio;

    @NotNull
    private final FittedRegionFactory fittedRegionFactory;
    private final ExecutorService executorService;
    private final Collection<SomaticVariant> variants;

    private final List<FittedPurity> all = Lists.newArrayList();
    private final List<ObservedRegion> filteredRegions = Lists.newArrayList();
    private final List<Double> ploidyRange;

    public FittedPurityFactory(final ExecutorService executorService, final CobaltChromosomes cobaltChromosomes, final double minPurity,
            final double maxPurity, final double purityIncrements, final double minPloidy, final double maxPloidy,
            final double somaticPenaltyWeight, final boolean tumorOnlyMode, @NotNull final FittedRegionFactory fittedRegionFactory,
            @NotNull final Collection<ObservedRegion> observedRegions, @NotNull final Collection<SomaticVariant> variants)
            throws ExecutionException, InterruptedException {
        this.executorService = executorService;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
        this.purityIncrements = purityIncrements;
        this.somaticPenaltyWeight = somaticPenaltyWeight;
        this.fittedRegionFactory = fittedRegionFactory;
        this.cobaltChromosomes = cobaltChromosomes;
        this.ploidyRange = ploidyRange(minPloidy, maxPloidy);

        final List<SomaticVariant> filteredVariants = Lists.newArrayList();
        final GenomePositionSelector<SomaticVariant> variantSelector = GenomePositionSelectorFactory.create(variants);

        int accumulatedBafCount = 0;
        double accumulatedWeightedRatio = 0;
        for (final ObservedRegion region : observedRegions) {
            if (useRegionToFitPurity(tumorOnlyMode, cobaltChromosomes, region)) {
                filteredRegions.add(region);
                variantSelector.select(region, filteredVariants::add);
                accumulatedBafCount += region.bafCount();
                accumulatedWeightedRatio += region.bafCount() * region.observedTumorRatio();
            }
        }

        this.totalBAFCount = accumulatedBafCount;
        this.averageFittingRatio = accumulatedWeightedRatio / accumulatedBafCount;
        this.variants = Downsample.downsample(MAX_SOMATICS_TO_FIT, filteredVariants);

        fitPurity();
    }

    @VisibleForTesting
    static boolean useRegionToFitPurity(boolean tumorOnlyMode, @NotNull final CobaltChromosomes cobaltChromosomes, @NotNull final ObservedRegion region) {
        if (region.bafCount() <= 0) {
            return false;
        }

        if (!positiveOrZero(region.observedTumorRatio())) {
            return false;
        }

        if (region.status() != GermlineStatus.DIPLOID) {
            return false;
        }

        if (Doubles.greaterThan(region.observedTumorRatio(), MAX_TUMOR_RATIO_TO_FIT)) {
            return false;
        }

        if (!cobaltChromosomes.contains(region.chromosome())) {
            return false;
        }

        CobaltChromosome chromosome = cobaltChromosomes.get(region.chromosome());
        if (tumorOnlyMode && chromosome.isAllosome()) {
            return false;
        }

        return chromosome.isNormal() && chromosome.isDiploid();
    }

    public List<FittedPurity> all() {
        return all;
    }

    private void fitPurity() throws ExecutionException, InterruptedException {
        final List<Future<List<FittedPurity>>> futures = Lists.newArrayList();
        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            futures.add(executorService.submit(callableFitPurity(purity)));
        }

        for (Future<List<FittedPurity>> future : futures) {
            List<FittedPurity> fittedPurities = future.get();

            if (!fittedPurities.isEmpty()) {
                all.addAll(fittedPurities);
            }
        }

        Collections.sort(all);
    }

    @NotNull
    private Callable<List<FittedPurity>> callableFitPurity(final double purity) {
        return () -> fitPurity(purity);
    }

    @NotNull
    private List<FittedPurity> fitPurity(final double purity) {
        final List<FittedPurity> fittedPurities = Lists.newArrayList();
        for (Double ploidy : ploidyRange) {
            double impliedNormFactor = PurityAdjuster.impliedNormFactor(averageFittingRatio, purity, ploidy);
            fittedPurities.add(fitPurity(purity, impliedNormFactor));
        }

        Collections.sort(fittedPurities);
        return fittedPurities;
    }

    private double weightWithBaf(double value, int bafCount) {
        return 1d * value * bafCount / totalBAFCount;
    }

    @NotNull
    private FittedPurity fitPurity(final double purity, final double normFactor) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double eventPenalty = 0;
        double deviationPenalty = 0;
        double diploidProportion = 0;
        double averagePloidy = 0;

        final List<FittedRegion> fittedRegions = Lists.newArrayList();
        for (final ObservedRegion enrichedRegion : filteredRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fitRegion(purity, normFactor, enrichedRegion);
            eventPenalty += weightWithBaf(fittedRegion.eventPenalty(), enrichedRegion.bafCount());
            deviationPenalty += weightWithBaf(fittedRegion.deviationPenalty(), enrichedRegion.bafCount());
            averagePloidy += weightWithBaf(fittedRegion.tumorCopyNumber(), enrichedRegion.bafCount());
            if (fittedRegion.isDiploid()) {
                diploidProportion += weightWithBaf(1, enrichedRegion.bafCount());
            }

            fittedRegions.add(fittedRegion);
        }

        final PurityAdjuster purityAdjuster = new PurityAdjusterAbnormalChromosome(purity, normFactor, cobaltChromosomes.chromosomes());
        final double somaticPenalty = Doubles.greaterThan(somaticPenaltyWeight, 0) ? somaticPenaltyWeight * SomaticPenaltyFactory.penalty(
                purityAdjuster,
                fittedRegions,
                variants) : 0;

        return builder.score(eventPenalty * deviationPenalty + somaticPenalty)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .somaticPenalty(somaticPenalty)
                .build();
    }

    @NotNull
    static List<Double> ploidyRange(double minPloidy, double maxPloidy) {
        List<Double> results = Lists.newArrayList();
        results.addAll(sequence(Math.max(0, minPloidy), Math.min(3, maxPloidy), 0.02));
        results.addAll(sequence(Math.max(3, minPloidy), Math.min(5, maxPloidy), 0.05));
        results.addAll(sequence(Math.max(5, minPloidy), maxPloidy, 0.1));
        results.add(maxPloidy);

        return results;
    }

    @NotNull
    private static List<Double> sequence(double inclusiveMin, double exclusiveMax, double increment) {
        List<Double> results = Lists.newArrayList();

        double ploidy = inclusiveMin;
        while (Doubles.lessThan(ploidy, exclusiveMax)) {
            results.add(ploidy);
            ploidy += increment;
        }
        return results;
    }

}
