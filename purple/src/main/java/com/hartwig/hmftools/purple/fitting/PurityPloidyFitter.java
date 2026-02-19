package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NORMAL;
import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleUtils.formatPurity;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.validateCopyNumbers;
import static com.hartwig.hmftools.purple.fitting.VariantPurityFitter.somaticFitIsWorse;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.purple.region.FittingRegion;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class PurityPloidyFitter
{
    private final SampleData mSampleData;
    private final RegionFitCalculator mRegionFitCalculator;
    private final List<ObservedRegion> mObservedRegions;
    private final Gender mGender;

    private final ExecutorService mExecutorService;
    private final PurpleConfig mConfig;
    private final boolean mTargetedMode;

    private final VariantPurityFitter mVariantPurityFitter;

    // output from fit
    private final List<PurpleCopyNumber> mCopyNumbers;
    private final List<ObservedRegion> mFittedRegions;

    private final List<FittedPurity> mCopyNumberFitCandidates;
    private FittedPurity mCopyNumberPurityFit;
    private FittedPurity mSomaticPurityFit;
    private FittedPurity mFinalPurityFit;
    private FittedPurityScore mFitPurityScore;
    private FittedPurityMethod mFitMethod;
    private BestFit mBestFit;
    private PurityAdjuster mPurityAdjuster;

    private final boolean mHasChimerism;
    private boolean mIsValid;

    public PurityPloidyFitter(
            final PurpleConfig config, final ReferenceData referenceData, final SampleData sampleData,
            final ExecutorService executorService,
            final RegionFitCalculator regionFitCalculator, final List<ObservedRegion> observedRegions, final Gender gender,
            final boolean hasChimerism)
    {
        mSampleData = sampleData;
        mConfig = config;
        mTargetedMode = referenceData.TargetRegions.hasTargetRegions();
        mExecutorService = executorService;
        mRegionFitCalculator = regionFitCalculator;
        mObservedRegions = observedRegions;
        mGender = gender;
        mHasChimerism = hasChimerism;

        mCopyNumbers = Lists.newArrayList();
        mFittedRegions = Lists.newArrayList();
        mCopyNumberFitCandidates = Lists.newArrayList();

        mVariantPurityFitter = new VariantPurityFitter(mConfig, referenceData, sampleData);
        mPurityAdjuster = null;

        mCopyNumberPurityFit = null;
        mFitPurityScore = null;
        mSomaticPurityFit = null;
        mFinalPurityFit = null;
        mBestFit = null;
        mFitMethod = null;
        mIsValid = true;
    }

    public List<PurpleCopyNumber> copyNumbers()
    {
        return mCopyNumbers;
    }

    public List<ObservedRegion> fittedRegions()
    {
        return mFittedRegions;
    }

    public boolean isValid()
    {
        return mIsValid;
    }

    public BestFit finalFit()
    {
        if(!mIsValid)
        {
            return null;
        }

        if(mBestFit == null)
        {
            mBestFit = new BestFit(mFinalPurityFit, mFitPurityScore, mFitMethod, mCopyNumberFitCandidates);
        }

        return mBestFit;
    }

    public PurityAdjuster purityAdjuster()
    {
        return mPurityAdjuster;
    }

    public void run()
    {
        mVariantPurityFitter.setState(mObservedRegions);

        performCopyNumberFit();

        if(mCopyNumberPurityFit == null)
        {
            PPL_LOGGER.error("failed to find copy number fit");
            mIsValid = false;
            return;
        }

        buildCopyNumbers(mCopyNumberPurityFit);

        performSomaticFit();

        if(mFinalPurityFit == null)
        {
            PPL_LOGGER.error("failed to find final fit");
            mIsValid = false;
            return;
        }

        if(mFinalPurityFit != mCopyNumberPurityFit)
        {
            buildCopyNumbers(mFinalPurityFit);
        }
    }

    private void performCopyNumberFit()
    {
        FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(
                mConfig, mExecutorService, mSampleData.Cobalt.CobaltChromosomes, mRegionFitCalculator, mObservedRegions,
                !mConfig.tumorOnlyMode() ? mVariantPurityFitter.fittingSomatics() : Collections.emptyList());

        if(!fittedPurityFactory.validDataForFit())
        {
            mIsValid = false;
            return;
        }

        try
        {
            fittedPurityFactory.fitPurity();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("error running copy number fit: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        mCopyNumberFitCandidates.addAll(fittedPurityFactory.getFittedPurities());

        Collections.sort(mCopyNumberFitCandidates);

        mCopyNumberPurityFit = mCopyNumberFitCandidates.get(0);

        List<FittedPurity> bestFitPerPurityCandidates = inRangeOfLowest(mCopyNumberPurityFit.score(), mCopyNumberFitCandidates);
        mFitPurityScore = FittedPurityScoreFactory.score(bestFitPerPurityCandidates);
    }

    private boolean noSignificantAneuploidy()
    {
        List<ObservedRegion> diploidRegions = mObservedRegions.stream()
                .filter(region -> useRegionToDetectAneuploidy(mSampleData.Cobalt.CobaltChromosomes, region))
                .toList();

        int highBafCount = 0;
        int totalBafCount = 0;
        for(ObservedRegion region : diploidRegions)
        {
            totalBafCount += region.bafCount();
            PPL_LOGGER.debug(format("AR: region %s:%d-%d has %d total points", region.chromosome(), region.start(), region.end(), region.bafCount()));
            if(region.observedBAF() > SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF
                    && region.bafCount() > SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT)
            {
                PPL_LOGGER.debug(format("AR: region contributes ratio: %d", region.bafCount()));
                highBafCount += region.bafCount();
            }
        }
        PPL_LOGGER.debug(format("AR: high diploid baf count: %d, total diploid baf count: %d", highBafCount, totalBafCount));
        double ratio = (double) highBafCount / totalBafCount;
        PPL_LOGGER.debug(format("aneuploidy ratio: %3f", ratio));
        return ratio < SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF;
    }

    private static boolean useRegionToDetectAneuploidy(final CobaltChromosomes cobaltChromosomes, final FittingRegion region)
    {
        if(region.bafCount() <= 0)
        {
            return false;
        }

        if(!positiveOrZero(region.observedTumorRatio()))
        {
            return false;
        }

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
        {
            return false;
        }

        return cobaltChromosomes.hasChromosome(region.chromosome());
    }

    private static boolean isCloseToDiploidAndHighPurity(final FittedPurity normalPurityFit)
    {
        return normalPurityFit.purity() > SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN
                && normalPurityFit.ploidy() > SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN
                && normalPurityFit.ploidy() < SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX;
    }

    private void performSomaticFit()
    {
        List<FittedPurity> diploidCandidates = BestFit.mostDiploidPerPurity(mCopyNumberFitCandidates);

        FittedPurity lowestPurityFit = !diploidCandidates.isEmpty() ?
                diploidCandidates.stream().min(Comparator.comparingDouble(FittedPurity::purity)).get() : mCopyNumberPurityFit;

        boolean highlyDiploid = isHighlyDiploid(mFitPurityScore); // max diploid proportion > 0.97
        if(mConfig.tumorOnlyMode() || mTargetedMode)
        {
            boolean noSignificantAneuploidy = noSignificantAneuploidy();
            boolean diploidHighPurity = isCloseToDiploidAndHighPurity(mCopyNumberPurityFit);
            if(mHasChimerism || (noSignificantAneuploidy && (diploidHighPurity || highlyDiploid)))
            {
                PPL_LOGGER.debug("Searching for somatic fit for non-chimeric diploid sample");
                mSomaticPurityFit = mVariantPurityFitter.calcSomaticOnlyFit(mCopyNumberFitCandidates);

                if(mSomaticPurityFit != null)
                {
                    mFinalPurityFit = mSomaticPurityFit;
                    mFitMethod = FittedPurityMethod.SOMATIC;
                    PPL_LOGGER.debug("Somatic fit found.");
                }
                else
                {
                    PPL_LOGGER.debug("Somatic fit not found, reverting to lowest purity fit.");
                    mFinalPurityFit = ImmutableFittedPurity.builder()
                            .purity(MIN_PURITY_DEFAULT)
                            .ploidy(2)
                            .normFactor(lowestPurityFit.normFactor())
                            .score(lowestPurityFit.score())
                            .diploidProportion(lowestPurityFit.diploidProportion())
                            .somaticPenalty(0) // defaults for the rest
                            .build();
                    mFitMethod = FittedPurityMethod.NO_TUMOR;
                }
            }
            else
            {
                mFinalPurityFit = mCopyNumberPurityFit;
                mFitMethod = FittedPurityMethod.NORMAL;
            }

            return;
        }

        boolean hasTumor = !highlyDiploid || mVariantPurityFitter.hasTumor();

        PPL_LOGGER.info("maxDiploidProportion({}) diploidCandidates({}) purityRange({} - {}) hasTumor({})",
                formatPurity(mFitPurityScore.maxDiploidProportion()), diploidCandidates.size(),
                formatPurity(mFitPurityScore.minPurity()), formatPurity(mFitPurityScore.maxPurity()), hasTumor);

        if(!hasTumor)
        {
            mFinalPurityFit = lowestPurityFit;
            mFitMethod = FittedPurityMethod.NO_TUMOR;
            return;
        }

        if(diploidCandidates.isEmpty())
        {
            PPL_LOGGER.warn("unable to use somatic fit as there are no diploid candidates");
            mFinalPurityFit = mCopyNumberPurityFit;
            mFitMethod = FittedPurityMethod.NORMAL;
            return;
        }

        boolean exceedsPuritySpread = Doubles.greaterOrEqual(mFitPurityScore.puritySpread(), mConfig.SomaticFitting.MinSomaticPuritySpread);

        boolean useSomatics = exceedsPuritySpread && highlyDiploid;

        if(!useSomatics)
        {
            mFinalPurityFit = mCopyNumberPurityFit;
            mFitMethod = FittedPurityMethod.NORMAL;
            return;
        }

        mSomaticPurityFit = mVariantPurityFitter.calcSomaticFit(diploidCandidates, mCopyNumbers, mGender);

        if(mSomaticPurityFit == null)
        {
            mFinalPurityFit = lowestPurityFit;
            mFitMethod = FittedPurityMethod.NO_TUMOR;
        }
        else
        {
            if(!somaticFitIsWorse(mCopyNumberPurityFit, mSomaticPurityFit, mConfig.SomaticFitting))
            {
                mFinalPurityFit = mSomaticPurityFit;
                mFitMethod = FittedPurityMethod.SOMATIC;
            }
            else
            {
                mFinalPurityFit = mCopyNumberPurityFit;
                mFitMethod = FittedPurityMethod.NORMAL;
            }
        }
    }

    private void buildCopyNumbers(final FittedPurity fittedPurity)
    {
        mCopyNumbers.clear();
        mFittedRegions.clear();

        mPurityAdjuster = new PurityAdjuster(fittedPurity.purity(), fittedPurity.normFactor(), mSampleData.Cobalt.CobaltChromosomes);

        final AmberData amberData = mSampleData.Amber;
        final CobaltData cobaltData = mSampleData.Cobalt;

        PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(
                mConfig.Fitting.MinDiploidTumorRatioCount,
                mConfig.Fitting.MinDiploidTumorRatioCountAtCentromere,
                amberData.AverageTumorDepth,
                fittedPurity.ploidy(),
                mPurityAdjuster,
                cobaltData.CobaltChromosomes);

        PPL_LOGGER.debug("building copy numbers");
        mFittedRegions.addAll(mRegionFitCalculator.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), mObservedRegions));

        copyNumberFactory.buildCopyNumbers(mFittedRegions, mSampleData.SvCache.somaticVariants());

        mCopyNumbers.addAll(copyNumberFactory.copyNumbers());

        if(!validateCopyNumbers(mCopyNumbers))
        {
            PPL_LOGGER.warn("invalid copy numbers");
            mIsValid = false;
        }
    }

    /*
    We trigger SOMAITC mode if
    [% of BAF sites >0.6] < 0.008
    && [
    MaxDiploidProportion > 0.97
    | (purity > 0.92 && ploidy =~2)]"
     */
    private boolean isHighlyDiploid(final FittedPurityScore score)
    {
        return Doubles.greaterOrEqual(score.maxDiploidProportion(), mConfig.SomaticFitting.HighlyDiploidPercentage);
    }

    private static List<FittedPurity> inRangeOfLowest(double lowestScore, final List<FittedPurity> purities)
    {
        return purities.stream().filter(inRangeOfLowest(lowestScore)).collect(toList());
    }

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;

    private static Predicate<FittedPurity> inRangeOfLowest(final double score)
    {
        return fittedPurity ->
        {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }

    public static BestFit buildGermlineBestFit()
    {
        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(1)
                .ploidy(2)
                .normFactor(1)
                .diploidProportion(1)
                .somaticPenalty(0)
                .score(0)
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPloidy(fittedPurity.ploidy())
                .maxPloidy(fittedPurity.ploidy())
                .minPurity(fittedPurity.purity())
                .maxPurity(fittedPurity.purity())
                .minDiploidProportion(1)
                .maxDiploidProportion(1)
                .build();

        return new BestFit(fittedPurity, score, NORMAL, List.of(fittedPurity));
    }
}
