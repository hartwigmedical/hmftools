package com.hartwig.hmftools.purple.fitting;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NORMAL;
import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_CONTAMINATION_CUTOFF;
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
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory;
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

    private final double ContaminationLevel;
    private boolean mIsValid;

    public PurityPloidyFitter(
            final PurpleConfig config, final ReferenceData referenceData, final SampleData sampleData,
            final ExecutorService executorService,
            final RegionFitCalculator regionFitCalculator, final List<ObservedRegion> observedRegions, final Gender gender)
    {
        mSampleData = sampleData;
        mConfig = config;
        mTargetedMode = referenceData.TargetRegions.hasTargetRegions();
        mExecutorService = executorService;
        mRegionFitCalculator = regionFitCalculator;
        mObservedRegions = observedRegions;
        mGender = gender;
        ContaminationLevel = sampleData.Amber.Contamination;

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
            AneuploidyDetector aneuploidyDetector =
                    new AneuploidyDetector(mObservedRegions, mSampleData.Amber.ChromosomeBafs, mSampleData.Cobalt.CobaltChromosomes);
            boolean noSignificantAneuploidy = !aneuploidyDetector.hasAneuploidy();
            //            boolean noSignificantAneuploidy = noSignificantAneuploidy();
            boolean diploidHighPurity = isCloseToDiploidAndHighPurity(mCopyNumberPurityFit);
            final boolean significantContamination = ContaminationLevel > SOMATIC_FIT_CONTAMINATION_CUTOFF;
            if(significantContamination || (noSignificantAneuploidy && (diploidHighPurity || highlyDiploid)))
            {
                PPL_LOGGER.debug("Searching for somatic fit for non-chimeric diploid sample. Significant contamination: {}", significantContamination);
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
                    mFinalPurityFit =
                            new FittedPurity(MIN_PURITY_DEFAULT, lowestPurityFit.normFactor(), 2, lowestPurityFit.score(), lowestPurityFit.diploidProportion(), 0);
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
        FittedPurity fittedPurity = new FittedPurity(1.0, 1.0, 2.0, 0.0, 1.0, 0.0);

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
