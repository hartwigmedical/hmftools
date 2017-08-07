package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory {

    private final int maxPloidy;
    private final double minPurity;
    private final double maxPurity;
    private final double purityIncrements;
    private final double minNormFactor;
    private final double maxNormFactor;
    private final double normFactorIncrements;
    @NotNull
    private final FittedRegionFactory fittedRegionFactory;
    private final ExecutorService executorService;

    private final List<FittedPurity> allPurities = Lists.newArrayList();
    private final List<FittedPurity> bestScoringPerPurity = Lists.newArrayList();

    public FittedPurityFactory(final ExecutorService executorService, final int maxPloidy, final double minPurity, final double maxPurity, final double purityIncrements,
            final double minNormFactor, final double maxNormFactor, final double normFactorIncrements,
            @NotNull final FittedRegionFactory fittedRegionFactory, @NotNull final Collection<ObservedRegion> observedRegions)
            throws ExecutionException, InterruptedException {
        this.executorService = executorService;
        this.maxPloidy = maxPloidy;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
        this.purityIncrements = purityIncrements;
        this.minNormFactor = minNormFactor;
        this.maxNormFactor = maxNormFactor;
        this.normFactorIncrements = normFactorIncrements;
        this.fittedRegionFactory = fittedRegionFactory;

        fitPurity(observedRegions);
    }

    public Optional<FittedPurity> bestFit() {
        return allPurities.isEmpty() ? Optional.empty() : Optional.of(allPurities.get(0));
    }

    public List<FittedPurity> allFits() {
        return allPurities;
    }

    public List<FittedPurity> bestFitPerPurity() {
        return bestScoringPerPurity;
    }

    private void fitPurity(@NotNull final Collection<ObservedRegion> observedRegions) throws ExecutionException, InterruptedException {
        int totalBAFCount = 0;
        final List<ObservedRegion> filteredRegions = Lists.newArrayList();

        for (final ObservedRegion region : observedRegions) {
            final Chromosome chromosome = HumanChromosome.valueOf(region);

            if (region.bafCount() > 0 && positiveOrZero(region.observedTumorRatio()) && chromosome.isAutosome()
                    && region.status() == FreecStatus.SOMATIC) {
                totalBAFCount += region.bafCount();
                filteredRegions.add(region);
            }
        }

        List<Future<List<FittedPurity>>> futures = Lists.newArrayList();
        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            futures.add(executorService.submit(callableFitPurity(purity, totalBAFCount, filteredRegions)));
        }

        for (Future<List<FittedPurity>> future : futures) {
            List<FittedPurity> fittedPurities = future.get();

            if (!fittedPurities.isEmpty()) {
                bestScoringPerPurity.add(fittedPurities.get(0));
            }
            allPurities.addAll(fittedPurities);
        }
        Collections.sort(bestScoringPerPurity);
        Collections.sort(allPurities);
    }

    private Callable<List<FittedPurity>> callableFitPurity(final double purity, final double totalBAFCount,
            @NotNull final Collection<ObservedRegion> observedRegions) {
        return () -> fitPurity(purity, totalBAFCount, observedRegions);
    }

    private List<FittedPurity> fitPurity(final double purity, final double totalBAFCount,
            @NotNull final Collection<ObservedRegion> observedRegions) {

        final List<FittedPurity> fittedPurities = Lists.newArrayList();
        for (double normFactor = minNormFactor; lessOrEqual(normFactor, maxNormFactor); normFactor += normFactorIncrements) {
            double impliedPloidy = PurityAdjuster.impliedSamplePloidy(purity, normFactor);

            if (greaterOrEqual(impliedPloidy, 1) && lessOrEqual(impliedPloidy, maxPloidy)) {
                fittedPurities.add(fitPurity(purity, normFactor, totalBAFCount, observedRegions));
            }
        }

        Collections.sort(fittedPurities);
        return fittedPurities;
    }

    @NotNull
    private FittedPurity fitPurity(final double purity, final double normFactor, final double totalBafCount,
            @NotNull final Collection<ObservedRegion> observedRegions) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0;
        double diploidProportion = 0;
        double modelBAFDeviation = 0;
        double averagePloidy = 0;

        for (final ObservedRegion enrichedRegion : observedRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fitRegion(purity, normFactor, enrichedRegion);
            modelDeviation += enrichedRegion.bafCount() / totalBafCount * fittedRegion.deviation();
            modelBAFDeviation += enrichedRegion.bafCount() / totalBafCount * fittedRegion.bafDeviation();
            averagePloidy += fittedRegion.tumorCopyNumber() * fittedRegion.bafCount() / totalBafCount;
            if (fittedRegion.fittedPloidy() == 2) {
                diploidProportion += enrichedRegion.bafCount() / totalBafCount;
            }
        }

        return builder.score(modelDeviation)
                .modelBAFDeviation(modelBAFDeviation)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .build();
    }
}
