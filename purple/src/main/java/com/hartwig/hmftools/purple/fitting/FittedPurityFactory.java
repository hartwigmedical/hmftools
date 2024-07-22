package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.purple.FittingConfig;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.fittingsnv.SomaticDeviation;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Downsample;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class FittedPurityFactory
{
    private final PurpleConfig mConfig;
    private final double mSomaticPenaltyWeight;
    private final Map<String,Double> mObservedRatioMap;

    private final int mTotalBAFCount;
    private final double mAverageFittingRatio;

    private final RegionFitCalculator mRegionFitCalculator;
    private final ExecutorService mExecutorService;

    private final List<FittedPurity> mFittedPurities;
    private final List<ObservedRegionData> mFilteredObservedRegions;
    private final List<Double> mPloidyRange;

    private static final int MAX_SOMATICS_TO_FIT = 1000;
    private static final double MAX_TUMOR_RATIO_TO_FIT = 3;

    public FittedPurityFactory(
            final PurpleConfig config, final ExecutorService executorService, final CobaltChromosomes cobaltChromosomes,
            final RegionFitCalculator regionFitCalculator, final Collection<ObservedRegion> observedRegions, final List<SomaticVariant> variants)
    {
        mConfig = config;
        mExecutorService = executorService;

        mSomaticPenaltyWeight = config.SomaticFitting.SomaticPenaltyWeight;

        mRegionFitCalculator = regionFitCalculator;
        mObservedRatioMap = cobaltChromosomes.chromosomes().stream()
                .collect(Collectors.toMap(CobaltChromosome::contig, CobaltChromosome::actualRatio));

        mPloidyRange = ploidyRange(mConfig.Fitting.MinPloidy, mConfig.Fitting.MaxPloidy);

        mFittedPurities = Lists.newArrayList();
        mFilteredObservedRegions = Lists.newArrayList();

        final List<SomaticVariant> filteredVariants = Lists.newArrayList();
        final GenomePositionSelector<SomaticVariant> variantSelector = GenomePositionSelectorFactory.create(variants);

        boolean tumorOnlyMode = mConfig.tumorOnlyMode();
        int accumulatedBafCount = 0;
        double accumulatedWeightedRatio = 0;

        for(ObservedRegion region : observedRegions)
        {
            if(useRegionToFitPurity(tumorOnlyMode, cobaltChromosomes, region))
            {
                variantSelector.select(region, filteredVariants::add);

                ObservedRegionData observedRegion = new ObservedRegionData(region);
                mFilteredObservedRegions.add(observedRegion);

                accumulatedBafCount += region.bafCount();
                accumulatedWeightedRatio += region.bafCount() * region.observedTumorRatio();
            }
        }

        mTotalBAFCount = accumulatedBafCount;
        mAverageFittingRatio = accumulatedWeightedRatio / accumulatedBafCount;

        List<SomaticVariant> downsampleVariants = Downsample.downsample(MAX_SOMATICS_TO_FIT, filteredVariants);

        // assign down-sampled variants to the filtered regions
        final GenomePositionSelector<SomaticVariant> filteredVariantSelector = GenomePositionSelectorFactory.create(downsampleVariants);

        for(ObservedRegionData regionData : mFilteredObservedRegions)
        {
            filteredVariantSelector.select(regionData.Region, regionData::addVariant);
        }
    }

    public List<FittedPurity> getFittedPurities() { return mFittedPurities; }

    public void fitPurity() throws ExecutionException, InterruptedException
    {
        FittingConfig config = mConfig.Fitting;

        if(mConfig.Threads <= 1)
        {
            for(double purity = config.MinPurity; lessOrEqual(purity, config.MaxPurity); purity += config.PurityIncrement)
            {
                mFittedPurities.addAll(fitPurity(purity));
            }
        }
        else
        {
            List<Future<List<FittedPurity>>> futures = Lists.newArrayList();
            for(double purity = config.MinPurity; lessOrEqual(purity, config.MaxPurity); purity += config.PurityIncrement)
            {
                futures.add(mExecutorService.submit(callableFitPurity(purity)));
            }

            for(Future<List<FittedPurity>> future : futures)
            {
                List<FittedPurity> fittedPurities = future.get();

                if(!fittedPurities.isEmpty())
                {
                    mFittedPurities.addAll(fittedPurities);
                }
            }
        }

        Collections.sort(mFittedPurities);
    }

    private Callable<List<FittedPurity>> callableFitPurity(final double purity)
    {
        return () -> fitPurity(purity);
    }

    private List<FittedPurity> fitPurity(final double purity)
    {
        final List<FittedPurity> fittedPurities = Lists.newArrayList();
        for(Double ploidy : mPloidyRange)
        {
            double impliedNormFactor = PurityAdjuster.impliedNormFactor(mAverageFittingRatio, purity, ploidy);
            fittedPurities.add(fitPurity(purity, impliedNormFactor));
        }

        return fittedPurities;
    }

    private FittedPurity fitPurity(final double purity, final double normFactor)
    {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double eventPenalty = 0;
        double deviationPenalty = 0;
        double diploidProportion = 0;
        double averagePloidy = 0;

        double somaticPenaltyTotal = 0;
        int somaticVariantCount = 0;
        final SomaticDeviation somaticDeviation = SomaticDeviation.INSTANCE;
        PurityAdjuster purityAdjuster = new PurityAdjuster(mObservedRatioMap, purity, normFactor);

        for(ObservedRegionData regionData : mFilteredObservedRegions)
        {
            ObservedRegion region = regionData.Region;
            RegionFitCalcs regionFitCalcs = mRegionFitCalculator.calculateRegionFit(purity, normFactor, region);

            int bafCount = region.bafCount();
            eventPenalty += weightWithBaf(regionFitCalcs.EventPenalty, bafCount);
            deviationPenalty += weightWithBaf(regionFitCalcs.DeviationPenalty, bafCount);
            averagePloidy += weightWithBaf(regionFitCalcs.TumorCopyNumber, bafCount);

            if(regionFitCalcs.isDiploid())
            {
                diploidProportion += weightWithBaf(1, bafCount);
            }

            for(SomaticVariant variant : regionData.Variants)
            {
                ++somaticVariantCount;

                double variantPenalty = somaticDeviation.deviationFromMax(
                        purityAdjuster, region.chromosome(), regionFitCalcs.majorAlleleCopyNumber(), regionFitCalcs.TumorCopyNumber, variant);

                somaticPenaltyTotal += mSomaticPenaltyWeight * variantPenalty;
            }

            /*
            PPL_LOGGER.trace(format("region(%s:%d-%d) fit(purity=%.2f norm=%.4f) somaticPenTotal(%.4f) devPen(%.4f) eventPen(%.4f) avgPloidy(%.4f)",
                    region.chromosome(), region.start(), region.end(),
                    purity, normFactor, somaticPenaltyTotal, deviationPenalty, eventPenalty, averagePloidy));
            */
        }

        double somaticPenalty = mSomaticPenaltyWeight > 0 && somaticVariantCount > 0 ? somaticPenaltyTotal / somaticVariantCount : 0;

        return builder.score(eventPenalty * deviationPenalty + somaticPenalty)
                .diploidProportion(diploidProportion)
                .ploidy(averagePloidy)
                .somaticPenalty(somaticPenalty)
                .build();
    }

    public static RegionFitCalculator createFittedRegionFactory(
            final int averageTumorDepth, final CobaltChromosomes cobaltChromosomes, final FittingConfig fitScoreConfig)
    {
        return new RegionFitCalculator(cobaltChromosomes, fitScoreConfig, averageTumorDepth);
    }

    private static boolean useRegionToFitPurity(boolean tumorOnlyMode, final CobaltChromosomes cobaltChromosomes, final ObservedRegion region)
    {
        if(region.bafCount() <= 0)
            return false;

        if(!positiveOrZero(region.observedTumorRatio()))
            return false;

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
            return false;

        if(Doubles.greaterThan(region.observedTumorRatio(), MAX_TUMOR_RATIO_TO_FIT))
            return false;

        if(!cobaltChromosomes.hasChromosome(region.chromosome()))
            return false;

        CobaltChromosome chromosome = cobaltChromosomes.get(region.chromosome());
        if(tumorOnlyMode && chromosome.isAllosome())
            return false;

        return chromosome.isNormal() && chromosome.isDiploid();
    }

    private double weightWithBaf(double value, int bafCount)
    {
        return 1d * value * bafCount / mTotalBAFCount;
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
