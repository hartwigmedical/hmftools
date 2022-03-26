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
import com.hartwig.hmftools.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.purple.somatic.SomaticData;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory
{
    private final double mMinPurity;
    private final double mMaxPurity;
    private final double mPurityIncrements;
    private final double mSomaticPenaltyWeight;
    private final CobaltChromosomes mCobaltChromosomes;

    private final int mTotalBAFCount;
    private final double mAverageFittingRatio;

    @NotNull
    private final FittedRegionFactory mFittedRegionFactory;
    private final ExecutorService mExecutorService;
    private final List<SomaticData> mVariants;

    private final List<FittedPurity> mAll = Lists.newArrayList();
    private final List<ObservedRegion> mFilteredRegions = Lists.newArrayList();
    private final List<Double> mPloidyRange;

    private static final int MAX_SOMATICS_TO_FIT = 1000;
    private static final double MAX_TUMOR_RATIO_TO_FIT = 3;

    public FittedPurityFactory(
            final ExecutorService executorService, final CobaltChromosomes cobaltChromosomes, final double minPurity,
            final double maxPurity, final double purityIncrements, final double minPloidy, final double maxPloidy,
            final double somaticPenaltyWeight, final boolean tumorOnlyMode, final FittedRegionFactory fittedRegionFactory,
            final Collection<ObservedRegion> observedRegions, final List<SomaticData> variants)
            throws ExecutionException, InterruptedException
    {
        mExecutorService = executorService;
        mMinPurity = minPurity;
        mMaxPurity = maxPurity;
        mPurityIncrements = purityIncrements;
        mSomaticPenaltyWeight = somaticPenaltyWeight;
        mFittedRegionFactory = fittedRegionFactory;
        mCobaltChromosomes = cobaltChromosomes;
        mPloidyRange = ploidyRange(minPloidy, maxPloidy);

        final List<SomaticData> filteredVariants = Lists.newArrayList();
        final GenomePositionSelector<SomaticData> variantSelector = GenomePositionSelectorFactory.create(variants);

        int accumulatedBafCount = 0;
        double accumulatedWeightedRatio = 0;
        for(final ObservedRegion region : observedRegions)
        {
            if(useRegionToFitPurity(tumorOnlyMode, cobaltChromosomes, region))
            {
                mFilteredRegions.add(region);
                variantSelector.select(region, filteredVariants::add);
                accumulatedBafCount += region.bafCount();
                accumulatedWeightedRatio += region.bafCount() * region.observedTumorRatio();
            }
        }

        mTotalBAFCount = accumulatedBafCount;
        mAverageFittingRatio = accumulatedWeightedRatio / accumulatedBafCount;
        mVariants = Downsample.downsample(MAX_SOMATICS_TO_FIT, filteredVariants);

        fitPurity();
    }

    @VisibleForTesting
    static boolean useRegionToFitPurity(boolean tumorOnlyMode, final CobaltChromosomes cobaltChromosomes, final ObservedRegion region)
    {
        if(region.bafCount() <= 0)
            return false;

        if(!positiveOrZero(region.observedTumorRatio()))
            return false;

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
            return false;

        if(Doubles.greaterThan(region.observedTumorRatio(), MAX_TUMOR_RATIO_TO_FIT))
            return false;

        if(!cobaltChromosomes.contains(region.chromosome()))
            return false;

        CobaltChromosome chromosome = cobaltChromosomes.get(region.chromosome());
        if(tumorOnlyMode && chromosome.isAllosome())
            return false;

        return chromosome.isNormal() && chromosome.isDiploid();
    }

    public List<FittedPurity> all()
    {
        return mAll;
    }

    private void fitPurity() throws ExecutionException, InterruptedException
    {
        final List<Future<List<FittedPurity>>> futures = Lists.newArrayList();
        for(double purity = mMinPurity; lessOrEqual(purity, mMaxPurity); purity += mPurityIncrements)
        {
            futures.add(mExecutorService.submit(callableFitPurity(purity)));
        }

        for(Future<List<FittedPurity>> future : futures)
        {
            List<FittedPurity> fittedPurities = future.get();

            if(!fittedPurities.isEmpty())
            {
                mAll.addAll(fittedPurities);
            }
        }

        Collections.sort(mAll);
    }

    @NotNull
    private Callable<List<FittedPurity>> callableFitPurity(final double purity)
    {
        return () -> fitPurity(purity);
    }

    @NotNull
    private List<FittedPurity> fitPurity(final double purity)
    {
        final List<FittedPurity> fittedPurities = Lists.newArrayList();
        for(Double ploidy : mPloidyRange)
        {
            double impliedNormFactor = PurityAdjuster.impliedNormFactor(mAverageFittingRatio, purity, ploidy);
            fittedPurities.add(fitPurity(purity, impliedNormFactor));
        }

        Collections.sort(fittedPurities);
        return fittedPurities;
    }

    private double weightWithBaf(double value, int bafCount)
    {
        return 1d * value * bafCount / mTotalBAFCount;
    }

    private FittedPurity fitPurity(final double purity, final double normFactor)
    {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double eventPenalty = 0;
        double deviationPenalty = 0;
        double diploidProportion = 0;
        double averagePloidy = 0;

        final List<FittedRegion> fittedRegions = Lists.newArrayList();
        for(final ObservedRegion enrichedRegion : mFilteredRegions)
        {
            final FittedRegion fittedRegion = mFittedRegionFactory.fitRegion(purity, normFactor, enrichedRegion);
            eventPenalty += weightWithBaf(fittedRegion.eventPenalty(), enrichedRegion.bafCount());
            deviationPenalty += weightWithBaf(fittedRegion.deviationPenalty(), enrichedRegion.bafCount());
            averagePloidy += weightWithBaf(fittedRegion.tumorCopyNumber(), enrichedRegion.bafCount());
            if(fittedRegion.isDiploid())
            {
                diploidProportion += weightWithBaf(1, enrichedRegion.bafCount());
            }

            fittedRegions.add(fittedRegion);
        }

        final PurityAdjuster purityAdjuster = new PurityAdjusterAbnormalChromosome(purity, normFactor, mCobaltChromosomes.chromosomes());

        final double somaticPenalty = Doubles.greaterThan(mSomaticPenaltyWeight, 0) ?
                mSomaticPenaltyWeight * SomaticPenaltyFactory.calcPenalty(purityAdjuster, fittedRegions, mVariants) : 0;

        return builder.score(eventPenalty * deviationPenalty + somaticPenalty)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .somaticPenalty(somaticPenalty)
                .build();
    }

    protected static List<Double> ploidyRange(double minPloidy, double maxPloidy)
    {
        List<Double> results = Lists.newArrayList();
        results.addAll(sequence(Math.max(0, minPloidy), Math.min(3, maxPloidy), 0.02));
        results.addAll(sequence(Math.max(3, minPloidy), Math.min(5, maxPloidy), 0.05));
        results.addAll(sequence(Math.max(5, minPloidy), maxPloidy, 0.1));
        results.add(maxPloidy);

        return results;
    }

    private static List<Double> sequence(double inclusiveMin, double exclusiveMax, double increment)
    {
        List<Double> results = Lists.newArrayList();

        double ploidy = inclusiveMin;
        while(Doubles.lessThan(ploidy, exclusiveMax))
        {
            results.add(ploidy);
            ploidy += increment;
        }
        return results;
    }

}
