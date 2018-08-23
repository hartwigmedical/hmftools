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
import com.hartwig.hmftools.common.collection.Downsample;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory {

    private static final int MAX_SOMATICS_TO_FIT = 1000;
    private static final double MAX_TUMOR_RATIO_TO_FIT = 3;

    private final int maxPloidy;
    private final Gender gender;
    private final double minPurity;
    private final double maxPurity;
    private final int totalBAFCount;
    private final double minNormFactor;
    private final double maxNormFactor;
    private final double purityIncrements;
    private final double normFactorIncrements;
    private final double somaticDeviationWeight;

    @NotNull
    private final FittedRegionFactory fittedRegionFactory;
    private final ExecutorService executorService;
    private final List<SomaticVariant> downsampledVariants;
    private final List<FittedPurity> bestScoringPerPurity = Lists.newArrayList();
    private final List<ObservedRegion> filteredObservedRegions = Lists.newArrayList();



    public FittedPurityFactory(final ExecutorService executorService, final Gender gender, final int maxPloidy, final double minPurity,
            final double maxPurity, final double purityIncrements, final double minNormFactor, final double maxNormFactor,
            final double normFactorIncrements, final double somaticDeviationWeight, @NotNull final FittedRegionFactory fittedRegionFactory,
            @NotNull final Collection<ObservedRegion> observedRegions, @NotNull final List<SomaticVariant> variants)
            throws ExecutionException, InterruptedException {
        this.executorService = executorService;
        this.maxPloidy = maxPloidy;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
        this.purityIncrements = purityIncrements;
        this.minNormFactor = minNormFactor;
        this.maxNormFactor = maxNormFactor;
        this.normFactorIncrements = normFactorIncrements;
        this.somaticDeviationWeight = somaticDeviationWeight;
        this.fittedRegionFactory = fittedRegionFactory;
        this.gender = gender;

        final List<SomaticVariant> filteredVariants = Lists.newArrayList();
        final GenomePositionSelector<SomaticVariant> variantSelector = GenomePositionSelectorFactory.create(variants);

        for (final ObservedRegion region : observedRegions) {
            final Chromosome chromosome = HumanChromosome.valueOf(region);
            if (region.bafCount() > 0 && positiveOrZero(region.observedTumorRatio()) && chromosome.isAutosome()
                    && region.status() == GermlineStatus.DIPLOID && Doubles.lessOrEqual(region.observedTumorRatio(),
                    MAX_TUMOR_RATIO_TO_FIT)) {
                filteredObservedRegions.add(region);
                variantSelector.select(region, filteredVariants::add);
            }
        }

        totalBAFCount = observedRegions.stream().mapToInt(ObservedRegion::bafCount).sum();
        downsampledVariants = Downsample.downsample(MAX_SOMATICS_TO_FIT, filteredVariants);

        fitPurity();
    }

    public List<FittedPurity> bestFitPerPurity() {
        return bestScoringPerPurity;
    }

    private void fitPurity() throws ExecutionException, InterruptedException {

        final List<Future<List<FittedPurity>>> futures = Lists.newArrayList();
        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            futures.add(executorService.submit(callableFitPurity(purity)));
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
    private Callable<List<FittedPurity>> callableFitPurity(final double purity) {
        return () -> fitPurity(purity);
    }

    @NotNull
    private List<FittedPurity> fitPurity(final double purity) {
        final List<FittedPurity> fittedPurities = Lists.newArrayList();
        for (double normFactor = minNormFactor; lessOrEqual(normFactor, maxNormFactor); normFactor += normFactorIncrements) {
            double impliedPloidy = PurityAdjuster.impliedSamplePloidy(purity, normFactor);

            if (greaterOrEqual(impliedPloidy, 1) && lessOrEqual(impliedPloidy, maxPloidy)) {
                fittedPurities.add(fitPurity(purity, normFactor));
            }
        }

        Collections.sort(fittedPurities);
        return fittedPurities;
    }

    @NotNull
    private FittedPurity fitPurity(final double purity, final double normFactor) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double ploidyPenalty = 0;
        double modelDeviation = 0;
        double diploidProportion = 0;
        double averagePloidy = 0;

        final List<FittedRegion> fittedRegions = Lists.newArrayList();
        for (final ObservedRegion enrichedRegion : filteredObservedRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fitRegion(purity, normFactor, enrichedRegion);
            ploidyPenalty += enrichedRegion.bafCount() / totalBAFCount * fittedRegion.ploidyPenalty();
            modelDeviation += enrichedRegion.bafCount() / totalBAFCount * fittedRegion.deviation();
            averagePloidy += fittedRegion.tumorCopyNumber() * fittedRegion.bafCount() / totalBAFCount;
            if (fittedRegion.modelPloidy() == 2) {
                diploidProportion += enrichedRegion.bafCount() / totalBAFCount;
            }

            fittedRegions.add(fittedRegion);
        }

        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);
        final double somaticDeviation = SomaticDeviationFactory.deviation(purityAdjuster, fittedRegions, downsampledVariants);

        return builder.score(ploidyPenalty * modelDeviation + somaticDeviationWeight * somaticDeviation)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .somaticDeviation(somaticDeviation)
                .build();
    }
}
