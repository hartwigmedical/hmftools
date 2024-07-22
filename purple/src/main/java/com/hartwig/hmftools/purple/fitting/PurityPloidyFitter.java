package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NORMAL;
import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleUtils.formatPurity;
import static com.hartwig.hmftools.purple.PurpleConstants.MAX_SOMATIC_FIT_DELETED_PERC;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.calculateDeletedDepthWindows;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.validateCopyNumbers;
import static com.hartwig.hmftools.purple.fitting.VariantPurityFitter.somaticFitIsWorse;
import static com.hartwig.hmftools.purple.fittingsnv.SomaticPurityFitter.useTumorOnlySomaticMode;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
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
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Segmentation;
import com.hartwig.hmftools.purple.sv.RecoverStructuralVariants;

public class PurityPloidyFitter
{
    private final SampleData mSampleData;
    private final RegionFitCalculator mRegionFitCalculator;
    private final List<ObservedRegion> mObservedRegions;

    private final ExecutorService mExecutorService;
    private final PurpleConfig mConfig;
    private final Segmentation mSegmentation;

    private VariantPurityFitter mVariantPurityFitter;

    // output from fit
    private final List<PurpleCopyNumber> mCopyNumbers;
    private final List<ObservedRegion> mFittedRegions;

    private List<FittedPurity> mCopyNumberFitCandidates;
    private FittedPurity mCopyNumberPurityFit;
    private FittedPurity mSomaticPurityFit;
    private FittedPurity mFinalPurityFit;
    private FittedPurityScore mFitPurityScore;
    private FittedPurityMethod mFitMethod;
    private BestFit mBestFit;
    private PurityAdjuster mPurityAdjuster;

    private boolean mIsValid;

    public PurityPloidyFitter(
            final PurpleConfig config, final ReferenceData referenceData, final SampleData sampleData, final ExecutorService executorService,
            final RegionFitCalculator regionFitCalculator, final List<ObservedRegion> observedRegions, final Segmentation segmentation)
    {
        mSampleData = sampleData;
        mConfig = config;
        mExecutorService = executorService;
        mRegionFitCalculator = regionFitCalculator;
        mObservedRegions = observedRegions;
        mSegmentation = segmentation;

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

    public List<PurpleCopyNumber> copyNumbers() { return mCopyNumbers; }
    public List<ObservedRegion> fittedRegions() { return mFittedRegions; }
    public List<FittedPurity> copyNumberFitCandidates() { return mCopyNumberFitCandidates; }
    public FittedPurity copyNumberFit() { return mCopyNumberPurityFit; }
    public FittedPurity somaticFit() { return mSomaticPurityFit; }
    public boolean isValid() { return mIsValid; }

    public BestFit finalFit()
    {
        if(!mIsValid)
            return null;

        if(mBestFit == null)
            mBestFit = new BestFit(mFinalPurityFit, mFitPurityScore, mFitMethod, mCopyNumberFitCandidates);

        return mBestFit;
    }

    public PurityAdjuster purityAdjuster() { return mPurityAdjuster; }

    public void run()
    {
        mVariantPurityFitter.setState(mObservedRegions);

        performCopyNumberFit();

        if(mCopyNumberPurityFit == null)
        {
            PPL_LOGGER.error("failed to find copy number fit");
            System.exit(1);
        }

        buildCopyNumbers(mCopyNumberPurityFit);

        performSomaticFit();

        if(mFinalPurityFit == null)
        {
            PPL_LOGGER.error("failed to find final fit");
            System.exit(1);
        }

        if(mFinalPurityFit != mCopyNumberPurityFit)
        {
            buildCopyNumbers(mFinalPurityFit);
        }

        determineFinalFit();
    }

    private void performCopyNumberFit()
    {
        FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(
                mConfig, mExecutorService, mSampleData.Cobalt.CobaltChromosomes, mRegionFitCalculator, mObservedRegions,
                mVariantPurityFitter.fittingSomatics());
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

        FittedPurity lowestScoreFit = mCopyNumberFitCandidates.get(0);

        mCopyNumberPurityFit = lowestScoreFit;

        List<FittedPurity> bestFitPerPurityCandidates = inRangeOfLowest(mCopyNumberPurityFit.score(), mCopyNumberFitCandidates);
        mFitPurityScore = FittedPurityScoreFactory.score(bestFitPerPurityCandidates);
    }

    private void performSomaticFit()
    {
        boolean exceedsPuritySpread = Doubles.greaterOrEqual(mFitPurityScore.puritySpread(), mConfig.SomaticFitting.MinSomaticPuritySpread);
        boolean highlyDiploid = isHighlyDiploid(mFitPurityScore);

        boolean hasTumor = !highlyDiploid || mVariantPurityFitter.hasTumor();
        List<FittedPurity> diploidCandidates = BestFit.mostDiploidPerPurity(mCopyNumberFitCandidates);

        PPL_LOGGER.info("maxDiploidProportion({}) diploidCandidates({}) purityRange({} - {}) hasTumor({})",
                formatPurity(mFitPurityScore.maxDiploidProportion()), diploidCandidates.size(),
                formatPurity(mFitPurityScore.minPurity()), formatPurity(mFitPurityScore.maxPurity()), hasTumor);

        FittedPurity lowestPurityFit = diploidCandidates.isEmpty() ?
                mCopyNumberPurityFit : diploidCandidates.stream().min(Comparator.comparingDouble(FittedPurity::purity)).get();

        // fit decision:
        // - if no tumor then take lowest score fit, method = NO_TUMOR, exit
        // - check for a tumor-only somatic fit

        if(!hasTumor)
        {
            mFinalPurityFit = mCopyNumberPurityFit;
            mFitMethod = FittedPurityMethod.NO_TUMOR;
            return;
        }

        if(mConfig.tumorOnlyMode() && useTumorOnlySomaticMode(mCopyNumberPurityFit))
        {
            mSomaticPurityFit = mVariantPurityFitter.tumorOnlySomaticFit(mCopyNumberFitCandidates);

            if(mSomaticPurityFit != null)
            {
                mFinalPurityFit = mSomaticPurityFit;
                mFitMethod = FittedPurityMethod.SOMATIC;
            }
            else
            {
                mFinalPurityFit = ImmutableFittedPurity.builder()
                        .purity(MIN_PURITY_DEFAULT).ploidy(2)
                        .score(0).diploidProportion(1).normFactor(1).somaticPenalty(0).build();
                mFitMethod = FittedPurityMethod.NO_TUMOR;
            }

            return;
        }

        if(diploidCandidates.isEmpty())
        {
            PPL_LOGGER.warn("unable to use somatic fit as there are no diploid candidates");
            mFinalPurityFit = mCopyNumberPurityFit;
            mFitMethod = FittedPurityMethod.NORMAL;
            return;
        }

        boolean useSomatics = mConfig.fitWithSomatics() && exceedsPuritySpread && highlyDiploid;

        if(!useSomatics)
        {
            mFinalPurityFit = mCopyNumberPurityFit;
            mFitMethod = FittedPurityMethod.NORMAL;
            return;
        }

        mSomaticPurityFit = mVariantPurityFitter.calcSomaticFit(diploidCandidates, mCopyNumbers);

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

    private void determineFinalFit()
    {
        if(mFitMethod != FittedPurityMethod.SOMATIC)
            return;

        // test the impact on deleted genes from a switch to use the somatic fit
        double deletedPercent = calculateDeletedDepthWindows(mCopyNumbers);

        if(deletedPercent >= MAX_SOMATIC_FIT_DELETED_PERC)
        {
            PPL_LOGGER.info(format("somatic fit(purity=%.3f ploidy=%.3f) deleted DW percent(%.3f), reverting to normal fit(purity=%.3f ploidy=%.3f)",
                    mSomaticPurityFit.purity(), mSomaticPurityFit.ploidy(), deletedPercent,
                    mCopyNumberPurityFit.purity(), mCopyNumberPurityFit.ploidy()));

            // re-build using the normal fit
            mFinalPurityFit = mCopyNumberPurityFit;

            buildCopyNumbers(mFinalPurityFit);
        }
        else
        {
            PPL_LOGGER.debug("somatic fit deleted DW percent({})", format("%.3f", deletedPercent));
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

        copyNumberFactory.buildCopyNumbers(mFittedRegions, mSampleData.SvCache.variants());

        int recoveredSVCount = RecoverStructuralVariants.recoverStructuralVariants(
                mSampleData, mConfig.SampleFiles, mConfig, mPurityAdjuster, copyNumberFactory.copyNumbers());

        if(recoveredSVCount > 0)
        {
            PPL_LOGGER.info("reapplying segmentation with {} recovered structural variants", recoveredSVCount);
            final List<ObservedRegion> recoveredObservedRegions =
                    mSegmentation.createObservedRegions(mSampleData.SvCache.variants(), amberData, cobaltData);

            PPL_LOGGER.info("recalculating copy number");
            mFittedRegions.clear();
            mFittedRegions.addAll(mRegionFitCalculator.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), recoveredObservedRegions));

            copyNumberFactory.buildCopyNumbers(mFittedRegions, mSampleData.SvCache.variants());
        }

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
