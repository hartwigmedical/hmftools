package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory {

    private static final double MAX_TUMOR_RATIO_TO_FIT = 3;

    private final int maxPloidy;
    private final double minPurity;
    private final double maxPurity;
    private final double purityIncrements;
    private final double minNormFactor;
    private final double maxNormFactor;
    private final double normFactorIncrements;
    private final Gender gender;
    @NotNull
    private final FittedRegionFactory fittedRegionFactory;
    private final ExecutorService executorService;
    private final Collection<SomaticVariant> variants;

    private final List<FittedPurity> bestScoringPerPurity = Lists.newArrayList();

    public FittedPurityFactory(final ExecutorService executorService, final Gender gender, final int maxPloidy, final double minPurity,
            final double maxPurity, final double purityIncrements, final double minNormFactor, final double maxNormFactor,
            final double normFactorIncrements, @NotNull final FittedRegionFactory fittedRegionFactory,
            @NotNull final Collection<ObservedRegion> observedRegions, @NotNull final Collection<SomaticVariant> variants)
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
        this.gender = gender;
        this.variants = variants;

        fitPurity(observedRegions);
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
                    && region.status() == GermlineStatus.DIPLOID && Doubles.lessOrEqual(region.observedTumorRatio(),
                    MAX_TUMOR_RATIO_TO_FIT)) {
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
        }

        Collections.sort(bestScoringPerPurity);
    }

    @NotNull
    private Callable<List<FittedPurity>> callableFitPurity(final double purity, final double totalBAFCount,
            @NotNull final Collection<ObservedRegion> observedRegions) {
        return () -> fitPurity(purity, totalBAFCount, observedRegions);
    }

    @NotNull
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
        double ploidyPenalty = 0;
        double modelDeviation = 0;
        double diploidProportion = 0;
        double averagePloidy = 0;

        //TODO: FIX AVERAGE PLOIDY

        final List<FittedRegion> fittedRegions = Lists.newArrayList();
        for (final ObservedRegion enrichedRegion : observedRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fitRegion(purity, normFactor, enrichedRegion);
            ploidyPenalty += enrichedRegion.bafCount() / totalBafCount * fittedRegion.ploidyPenalty();
            modelDeviation += enrichedRegion.bafCount() / totalBafCount * fittedRegion.deviation();
            averagePloidy += fittedRegion.tumorCopyNumber() * fittedRegion.bafCount() / totalBafCount;
            if (fittedRegion.modelPloidy() == 2) {
                diploidProportion += enrichedRegion.bafCount() / totalBafCount;
            }

            fittedRegions.add(fittedRegion);
        }

//        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);
//        double somaticDeviation = new SomaticDeviationFactory(purityAdjuster).deviation(fittedRegions, variants);


        return builder.score(ploidyPenalty * modelDeviation)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .somaticDeviation(0)
                .build();
    }
}
