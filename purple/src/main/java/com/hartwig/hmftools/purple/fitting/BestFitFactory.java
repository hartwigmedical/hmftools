package com.hartwig.hmftools.purple.fitting;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.purple.region.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleCommon.formatPurity;
import static com.hartwig.hmftools.purple.config.PurpleConstants.NO_TUMOR_BAF_TOTAL;
import static com.hartwig.hmftools.purple.config.PurpleConstants.NO_TUMOR_DEPTH_RATIO_MAX;
import static com.hartwig.hmftools.purple.config.PurpleConstants.NO_TUMOR_DEPTH_RATIO_MIN;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.purity.FittedPurityScoreFactory;
import com.hartwig.hmftools.common.purple.purity.ImmutableBestFit;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.purple.config.PurpleConfig;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class BestFitFactory
{
    private final PurpleConfig mConfig;
    private final SomaticPurityFitter mSomaticPurityFitter;

    private final int mMinReadCount;
    private final int mMaxReadCount;

    private int mSvHotspotCount;
    private int mSvFragmentReadCount;
    private int mSomaticHotspotCount;
    private int mAlleleReadCountTotal;
    private final List<SomaticVariant> mVariantsInReadCountRange;

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;

    @NotNull
    private final BestFit mBestFit;

    public BestFitFactory(
            final PurpleConfig config, int minReadCount, int maxReadCount,
            @NotNull final List<FittedPurity> allCandidates, final List<SomaticVariant> somatics,
            final List<StructuralVariant> structuralVariants, final List<ObservedRegion> observedRegions)
    {
        mConfig = config;
        mMinReadCount = minReadCount;
        mMaxReadCount = maxReadCount;

        mSvHotspotCount = 0;
        mSvFragmentReadCount = 0;
        mSomaticHotspotCount = 0;
        mAlleleReadCountTotal = 0;
        mVariantsInReadCountRange = Lists.newArrayList();

        mSomaticPurityFitter = new SomaticPurityFitter(
                config.SomaticFitting.MinPeakVariants, config.SomaticFitting.MinTotalVariants,
                config.Fitting.MinPurity, config.Fitting.MaxPurity);

        mBestFit = determineBestFit(allCandidates, somatics, structuralVariants, observedRegions);
    }

    @NotNull
    public BestFit bestFit() { return mBestFit; }

    private BestFit determineBestFit(
            final List<FittedPurity> allCandidates, final List<SomaticVariant> somatics,
            final List<StructuralVariant> structuralVariants, final List<ObservedRegion> observedRegions)
    {
        Collections.sort(allCandidates);
        FittedPurity lowestScoreFit = allCandidates.get(0);

        final List<FittedPurity> bestFitPerPurityCandidates = inRangeOfLowest(lowestScoreFit.score(), allCandidates);
        final FittedPurityScore score = FittedPurityScoreFactory.score(bestFitPerPurityCandidates);
        final ImmutableBestFit.Builder builder = ImmutableBestFit.builder().score(score).allFits(allCandidates);

        boolean exceedsPuritySpread = Doubles.greaterOrEqual(score.puritySpread(), mConfig.SomaticFitting.MinSomaticPuritySpread);
        boolean highlyDiploid = isHighlyDiploid(score);

        boolean hasTumor = !highlyDiploid || hasTumor(somatics, structuralVariants, observedRegions);
        final List<FittedPurity> diploidCandidates = BestFit.mostDiploidPerPurity(allCandidates);

        PPL_LOGGER.info("Sample maxDiploidProportion({}) diploidCandidates({}) purityRange({} - {}) hasTumor({})",
                formatPurity(score.maxDiploidProportion()), diploidCandidates.size(),
                formatPurity(score.minPurity()), formatPurity(score.maxPurity()), hasTumor);

        FittedPurity lowestPurityFit = diploidCandidates.isEmpty() ?
                lowestScoreFit : diploidCandidates.stream().min(Comparator.comparingDouble(FittedPurity::purity)).get();

        if(!hasTumor)
            return builder.fit(lowestPurityFit).method(FittedPurityMethod.NO_TUMOR).build();

        if(diploidCandidates.isEmpty())
        {
            PPL_LOGGER.warn("Unable to use somatic fit as there are no diploid candidates");
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();
        }

        boolean useSomatics = !mConfig.TumorOnlyMode && exceedsPuritySpread && highlyDiploid;

        if(!useSomatics)
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();

        final Optional<FittedPurity> somaticFit = mSomaticPurityFitter.fromSomatics(
                mVariantsInReadCountRange, structuralVariants, diploidCandidates);

        if(!somaticFit.isPresent())
        {
            return builder.fit(lowestPurityFit).method(FittedPurityMethod.SOMATIC).build();
        }
        else if(somaticFitIsWorse(lowestScoreFit, somaticFit.get()))
        {
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();
        }
        else
        {
            return builder.fit(somaticFit.get()).method(FittedPurityMethod.SOMATIC).build();
        }
    }

    private boolean hasTumor(
            final List<SomaticVariant> somatics, final List<StructuralVariant> structuralVariants, final List<ObservedRegion> observedRegions)
    {
        setSvSummary(structuralVariants);
        setSomaticSummary(somatics);

        if(mSomaticHotspotCount > 0 || mAlleleReadCountTotal >= mConfig.SomaticFitting.minTotalSomaticVariantAlleleReadCount())
        {
            PPL_LOGGER.info("Tumor evidence: somaticHotspotCount({}) alleleReadCountTotal({})",
                    mSomaticHotspotCount, mAlleleReadCountTotal);
            return true;
        }

        if(mSvHotspotCount > 0 || mSvFragmentReadCount >= mConfig.SomaticFitting.minTotalSvFragmentCount())
        {
            PPL_LOGGER.info("Tumor evidence: svHotspotCount({}) svFragmentReadCount({})",
                    mSvHotspotCount, mSvFragmentReadCount);
            return true;
        }

        int tumorEvidenceBafCountTotal = observedRegions.stream()
                .filter(x -> x.germlineStatus() == DIPLOID)
                .filter(x -> x.observedTumorRatio() < NO_TUMOR_DEPTH_RATIO_MIN || x.observedTumorRatio() > NO_TUMOR_DEPTH_RATIO_MAX)
                .mapToInt(x -> x.bafCount())
                .sum();

        if(tumorEvidenceBafCountTotal >= NO_TUMOR_BAF_TOTAL)
        {
            PPL_LOGGER.info("Tumor evidence: tumorEvidenceBafCountTotal({})", tumorEvidenceBafCountTotal);
            return true;
        }

        return false;
    }

    private boolean somaticFitIsWorse(final FittedPurity lowestScore, final FittedPurity somaticFit)
    {
        double lowestPurity = lowestScore.purity();
        double somaticPurity = somaticFit.purity();

        return Doubles.lessThan(lowestPurity, mConfig.SomaticFitting.MinSomaticPurity)
            && Doubles.lessThan(somaticPurity, mConfig.SomaticFitting.MinSomaticPurity) && Doubles.greaterThan(
                somaticPurity,
                lowestPurity);
    }

    private boolean isHighlyDiploid(@NotNull final FittedPurityScore score)
    {
        return Doubles.greaterOrEqual(score.maxDiploidProportion(), mConfig.SomaticFitting.HighlyDiploidPercentage);
    }

    private static List<FittedPurity> inRangeOfLowest(double lowestScore, @NotNull final List<FittedPurity> purities)
    {
        return purities.stream().filter(inRangeOfLowest(lowestScore)).collect(toList());
    }

    private static Predicate<FittedPurity> inRangeOfLowest(final double score)
    {
        return fittedPurity ->
        {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }

    private void setSvSummary(@NotNull final List<StructuralVariant> variants)
    {
        for(StructuralVariant variant : variants)
        {
            if(variant.isFiltered())
                continue;

            if(variant.hotspot())
                mSvHotspotCount++;

            Integer startTumorVariantFragmentCount = variant.start().tumorVariantFragmentCount();
            if(variant.end() != null && startTumorVariantFragmentCount != null)
            {
                mSvFragmentReadCount += startTumorVariantFragmentCount;
            }
        }
    }

    private void setSomaticSummary(final List<SomaticVariant> somatics)
    {
        for(SomaticVariant somatic : somatics)
        {
            if(somatic.isFiltered() || !somatic.isSnp())
                continue;

            if(somatic.isHotspot())
                mSomaticHotspotCount++;

            mAlleleReadCountTotal += somatic.alleleReadCount();

            if(somatic.totalReadCount() >= mMinReadCount && somatic.totalReadCount() <= mMaxReadCount)
            {
                mVariantsInReadCountRange.add(somatic);
            }
        }
    }
}
