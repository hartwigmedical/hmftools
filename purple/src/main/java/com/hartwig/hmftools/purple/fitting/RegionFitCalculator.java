package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.signum;

import static com.hartwig.hmftools.purple.PurpleConstants.BAF_PNT_5;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.purple.FittingConfig;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.region.PloidyDeviation;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.utils.Doubles;

public class RegionFitCalculator
{
    private final double mAmbiguousBaf;
    private final FittingConfig mFitScoreConfig;
    private final PloidyDeviation mPloidyDeviation;
    private final CobaltChromosomes mCobaltChromosomes;

    public RegionFitCalculator(final CobaltChromosomes cobaltChromosomes, final FittingConfig fitScoreConfig, int averageReadDepth)
    {
        mCobaltChromosomes = cobaltChromosomes;
        mFitScoreConfig = fitScoreConfig;

        mPloidyDeviation = new PloidyDeviation(
                fitScoreConfig.PloidyPenaltyStandardDeviation,
                fitScoreConfig.PloidyPenaltyMinStandardDeviationPerPloidy,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneMultiplier,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneAdditional,
                fitScoreConfig.PloidyPenaltyMinDeviation);

        mAmbiguousBaf = ExpectedBAF.expectedBAF(averageReadDepth);
    }

    public List<ObservedRegion> fitRegion(double purity, double normFactor, final Collection<ObservedRegion> observedRegions)
    {
        final Predicate<ObservedRegion> valid = observedRegion -> isAllowedRegion(mCobaltChromosomes, observedRegion);
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @VisibleForTesting
    public static boolean isAllowedRegion(final CobaltChromosomes cobaltChromosomes, final GenomeRegion region)
    {
        return cobaltChromosomes.hasChromosome(region.chromosome());
    }

    public ObservedRegion fitRegion(final double purity, final double normFactor, final ObservedRegion observedRegion)
    {
        RegionFitCalcs regionFitCalcs = calculateRegionFit(purity, normFactor, observedRegion);

        ObservedRegion fittedRegion = ObservedRegion.from(observedRegion);

        fittedRegion.setTumorCopyNumber(regionFitCalcs.TumorCopyNumber);
        fittedRegion.setTumorBAF(regionFitCalcs.TumorBAF);
        fittedRegion.setRefNormalisedCopyNumber(regionFitCalcs.RefNormalisedCopyNumber);
        fittedRegion.setMinorAlleleCopyNumberDeviation(regionFitCalcs.MinorAlleleCopyNumberDeviation);
        fittedRegion.setMajorAlleleCopyNumberDeviation(regionFitCalcs.MajorAlleleCopyNumberDeviation);
        fittedRegion.setDeviationPenalty(regionFitCalcs.DeviationPenalty);
        fittedRegion.setEventPenalty(regionFitCalcs.EventPenalty);

        return fittedRegion;
    }

    public RegionFitCalcs calculateRegionFit(final double purity, final double normFactor, final ObservedRegion observedRegion)
    {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purity, normFactor, mCobaltChromosomes);

        double observedTumorRatio = observedRegion.observedTumorRatio();

        double impliedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double observedBAF = observedRegion.observedBAF();
        double impliedBAF = impliedBaf(purityAdjuster, observedRegion.chromosome(), impliedCopyNumber, observedBAF);

        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());

        double majorAllelePloidy = impliedBAF * impliedCopyNumber;
        double minorAllelePloidy = impliedCopyNumber - majorAllelePloidy;

        double majorAllelePloidyDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAllelePloidy);
        double minorAllelePloidyDeviation = mPloidyDeviation.minorAlleleDeviation(purity, normFactor, minorAllelePloidy);

        double eventPenalty = calculateEventPenalty(mFitScoreConfig.PloidyPenaltyFactor, majorAllelePloidy, minorAllelePloidy);

        double deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF;

        if(mFitScoreConfig.GcRatioExponent > 0 || mFitScoreConfig.DeviationPenaltyGcMinAdjust > 0)
        {
            // NEW FORMULA:
            // deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF
            //  / max(DeviationPenaltyGcMinAdjust, observedTumorRatio^GcRatioExponent)

            double adjTumorRatio = mFitScoreConfig.GcRatioExponent > 0 ? pow(observedTumorRatio, mFitScoreConfig.GcRatioExponent) : 1;

            double deviationPenaltyDenom = max(mFitScoreConfig.DeviationPenaltyGcMinAdjust, adjTumorRatio);

            deviationPenalty /= deviationPenaltyDenom;
        }

        return new RegionFitCalcs(
                impliedCopyNumber, impliedBAF, Doubles.replaceNaNWithZero(refNormalisedCopyNumber),
                minorAllelePloidyDeviation, majorAllelePloidyDeviation, eventPenalty, deviationPenalty);
    }

    public static double calculateEventPenalty(double eventPenaltyFactor, double majorAllele, double minorAllele)
    {
        double wholeGenomeDoublingDistance = 1 + (abs(majorAllele - 2)) + (abs(minorAllele - 2));
        double singleEventDistance = (abs(majorAllele - 1)) + (abs(minorAllele - 1));

        return 1 + eventPenaltyFactor * min(singleEventDistance, wholeGenomeDoublingDistance);
    }

    private static final double MIN_CN_THRESHOLD = 0.1;

    private double impliedBaf(final PurityAdjuster purityAdjuster, final String chromosome, double copyNumber, double observedBAF)
    {
        if(!mCobaltChromosomes.hasChromosome(chromosome))
            return 1;

        CobaltChromosome cobaltChromosome = mCobaltChromosomes.get(chromosome);

        if(!cobaltChromosome.isNormal() || !cobaltChromosome.isDiploid() || Doubles.lessOrEqual(copyNumber, MIN_CN_THRESHOLD))
            return 1;

        if(Doubles.lessOrEqual(observedBAF, mAmbiguousBaf))
            return bafToMinimiseDeviation(purityAdjuster, chromosome, copyNumber, observedBAF);
        else
            return purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF);
    }

    @VisibleForTesting
    public double bafToMinimiseDeviation(final PurityAdjuster purityAdjuster, final String chromosome, double copyNumber, double observedBAF)
    {
        double minBAF = max(0, min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, BAF_PNT_5)));
        double maxBAF = max(0, min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF)));

        double estimatedBaf = estimateMinMaxBaf(copyNumber, minBAF, maxBAF);

        if(estimatedBaf != NO_CALC_BAF)
            return estimatedBaf;

        double majorAcnMin = minBAF * copyNumber;
        double majorAcnMax = maxBAF * copyNumber;
        double minorAcnMin = copyNumber - majorAcnMin;
        double minorAcnMax = copyNumber - majorAcnMax;

        double purity = purityAdjuster.purity();
        double normFactor = purityAdjuster.normFactor();

        // minimise
        double minBAFTotalDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAcnMin)
                        + mPloidyDeviation.minorAlleleDeviation(purity, normFactor, minorAcnMin);

        double maxBAFTotalDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAcnMax)
                        + mPloidyDeviation.minorAlleleDeviation(purity,normFactor, minorAcnMax);

        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? BAF_PNT_5 : observedBAF;
    }

    private static final double NO_CALC_BAF = -1;

    public static double estimateMinMaxBaf(double copyNumber, double minBAF, double maxBAF)
    {
        // calculate major and minor copy numbers from these min/max BAFs
        double majorAcnMin = minBAF * copyNumber;
        double majorAcnMax = maxBAF * copyNumber;
        double minorAcnMin = copyNumber - majorAcnMin;
        double minorAcnMax = copyNumber - majorAcnMax;

        // test for whole number solutions
        double minorAcnMaxCeil = ceil(minorAcnMax);

        boolean minorDiffIntegers = !Doubles.equal(
                signum(minorAcnMaxCeil - minorAcnMin), signum(minorAcnMaxCeil - minorAcnMax));

        double majorAcnMinCeil = ceil(majorAcnMin);

        boolean majorDiffIntegers = !Doubles.equal(
                signum(majorAcnMinCeil - majorAcnMin), signum(majorAcnMinCeil - majorAcnMax));

        if(!minorDiffIntegers && !majorDiffIntegers)
            return NO_CALC_BAF;

        // test for use of either only major or minor
        if(!minorDiffIntegers)
        {
            double bafFromMajor = majorAcnMinCeil / copyNumber;
            return bafFromMajor;
        }
        else if(!majorDiffIntegers)
        {
            double bafFromMinor = 1 - minorAcnMaxCeil / copyNumber;
            return bafFromMinor;
        }

        double halfCopyNumber = copyNumber * 0.5;
        double majorAcnMaxFloor = floor(majorAcnMax);
        double majorAcnEstimateHigh = majorAcnMaxFloor >= halfCopyNumber ? majorAcnMaxFloor : majorAcnMax;

        // test values for the major allele CN between the bounds of ..
        double majorAcnEstimateLow = majorAcnEstimateHigh - 1;

        double majorLowerBound = max(halfCopyNumber, halfCopyNumber);

        for(int minorInt = 0; minorInt < copyNumber; ++minorInt)
        {
            double majorAcn = copyNumber - minorInt;

            if(majorAcn > minorInt && majorAcn - minorInt <= 0.5)
            {
                majorAcnEstimateLow = majorAcn;
                break;
            }

            if(majorAcn >= majorLowerBound)
            {
                majorAcnEstimateLow = majorAcn;
            }
            else if(majorAcn < minorInt)
            {
                double previousFloor = floor(majorAcn + 1);

                if(previousFloor >= majorLowerBound)
                    majorAcnEstimateLow = previousFloor;

                break;
            }
        }

        majorAcnEstimateLow = min(majorAcnEstimateLow, majorAcnEstimateHigh);

        double minorAcnEstimateHigh = copyNumber - majorAcnEstimateHigh;
        double minorAcnEstimateLow = copyNumber - majorAcnEstimateLow;

        double diffHigh = majorAcnEstimateHigh - minorAcnEstimateHigh;
        double diffLow = majorAcnEstimateLow - minorAcnEstimateLow;

        double majorAcnEstimate = diffHigh <= diffLow ? majorAcnEstimateHigh : majorAcnEstimateLow;

        return max(majorAcnEstimate / copyNumber, BAF_PNT_5);
    }
}
