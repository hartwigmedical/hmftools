package com.hartwig.hmftools.purple.purity;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.signum;

import static com.hartwig.hmftools.purple.config.PurpleConstants.BAF_PNT_5;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.purple.config.FittingConfig;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.region.PloidyDeviation;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.utils.Doubles;

public class RegionFitCalculator
{
    private final double mAmbiguousBaf;
    private final double mPloidyPenaltyFactor;
    private final PloidyDeviation mPloidyDeviation;
    private final CobaltChromosomes mCobaltChromosomes;

    public RegionFitCalculator(final CobaltChromosomes cobaltChromosomes, final FittingConfig fitScoreConfig, int averageReadDepth)
    {
        mCobaltChromosomes = cobaltChromosomes;
        mPloidyPenaltyFactor = fitScoreConfig.PloidyPenaltyFactor;

        mPloidyDeviation = new PloidyDeviation(
                fitScoreConfig.PloidyPenaltyStandardDeviation,
                fitScoreConfig.PloidyPenaltyMinStandardDeviationPerPloidy,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneMultiplier,
                fitScoreConfig.PloidyPenaltyMajorAlleleSubOneAdditional,
                fitScoreConfig.PloidyPenaltyBaselineDeviation);

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

        double eventPenalty = calculateEventPenalty(mPloidyPenaltyFactor, majorAllelePloidy, minorAllelePloidy);
        double deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF;

        return new RegionFitCalcs(
                impliedCopyNumber, impliedBAF, Doubles.replaceNaNWithZero(refNormalisedCopyNumber),
                minorAllelePloidyDeviation, majorAllelePloidyDeviation, eventPenalty, deviationPenalty);
    }

    public static double calculateEventPenalty(double eventPenaltyFactor, double majorAllele, double minorAllele)
    {
        double wholeGenomeDoublingDistance = 1 + (abs(majorAllele - 2)) + (abs(minorAllele - 2));
        double singleEventDistance = (abs(majorAllele - 1)) + (abs(minorAllele - 1));

        return 1 + eventPenaltyFactor * Math.min(singleEventDistance, wholeGenomeDoublingDistance);
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

        /*
        return Doubles.lessOrEqual(observedBAF, mAmbiguousBaf)
                ? bafToMinimiseDeviation(purityAdjuster, chromosome, copyNumber, observedBAF)
                : purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF);
        */
    }

    @VisibleForTesting
    public double bafToMinimiseDeviation(final PurityAdjuster purityAdjuster, final String chromosome, double copyNumber, double observedBAF)
    {
        double minBAF = max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, BAF_PNT_5)));
        double maxBAF = max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF)));

        // calculate major and minor min/max copy number
        double majorAlleleCnMin = minBAF * copyNumber;
        double majorAlleleCnMax = maxBAF * copyNumber;
        double minorAlleleCnMin = copyNumber - majorAlleleCnMin;
        double minorAlleleCnMax = copyNumber - majorAlleleCnMax;

        // test for whole number solutions
        double minorAlleleCnMaxCeil = ceil(minorAlleleCnMax);

        boolean minorDiffIntegers = !Doubles.equal(
                signum(minorAlleleCnMaxCeil - minorAlleleCnMin), signum(minorAlleleCnMaxCeil - minorAlleleCnMax));

        double majorAlleleCnMinCeil = ceil(majorAlleleCnMin);

        boolean majorDiffIntegers = !Doubles.equal(
                signum(majorAlleleCnMinCeil - majorAlleleCnMin), signum(majorAlleleCnMinCeil - majorAlleleCnMax));

        if(minorDiffIntegers && majorDiffIntegers)
        {
            // select the solution which minimum difference between the major and minor CNs
            if(abs(majorAlleleCnMin - minorAlleleCnMin) < abs(majorAlleleCnMax - minorAlleleCnMax))
                return max(floor(majorAlleleCnMin) / copyNumber, BAF_PNT_5);
            else
                return majorAlleleCnMinCeil / copyNumber;
        }
        else if(majorDiffIntegers)
        {
            return majorAlleleCnMinCeil / copyNumber;
        }
        else if(minorDiffIntegers)
        {
            return 1 - minorAlleleCnMaxCeil / copyNumber;
        }

        double purity = purityAdjuster.purity();
        double normFactor = purityAdjuster.normFactor();

        // minimise
        double minBAFTotalDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAlleleCnMin)
                        + mPloidyDeviation.minorAlleleDeviation(purity, normFactor, minorAlleleCnMin);

        double maxBAFTotalDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAlleleCnMax)
                        + mPloidyDeviation.minorAlleleDeviation(purity,normFactor, minorAlleleCnMax);

        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? BAF_PNT_5 : observedBAF;
    }

    private double bafToMinimiseDeviationOld(final PurityAdjuster purityAdjuster, final String chromosome, double impliedCopyNumber)
    {
        final double minBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, 0.5)));
        final double maxBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, mAmbiguousBaf)));

        // Major Ploidy
        final double minBAFMajorAllelePloidy = minBAF * impliedCopyNumber;
        final double maxBAFMajorAllelePloidy = maxBAF * impliedCopyNumber;

        // Major Ploidy crosses whole number?
        final double minBAFMajorAllelePloidyCeil = Math.ceil(minBAFMajorAllelePloidy);
        if(!Doubles.equal(Math.signum(minBAFMajorAllelePloidyCeil - minBAFMajorAllelePloidy),
                Math.signum(minBAFMajorAllelePloidyCeil - maxBAFMajorAllelePloidy)))
        {
            return minBAFMajorAllelePloidyCeil / impliedCopyNumber;
        }

        // Minor Ploidy
        final double minBAFMinorAllelePloidy = impliedCopyNumber - minBAFMajorAllelePloidy;
        final double maxBAFMinorAllelePloidy = impliedCopyNumber - maxBAFMajorAllelePloidy;

        // Minor Ploidy crosses whole number?
        final double maxBAFMinorAllelePloidyCeil = Math.ceil(maxBAFMinorAllelePloidy);
        if(!Doubles.equal(Math.signum(maxBAFMinorAllelePloidyCeil - minBAFMinorAllelePloidy),
                Math.signum(maxBAFMinorAllelePloidyCeil - maxBAFMinorAllelePloidy)))
        {
            return 1 - maxBAFMinorAllelePloidyCeil / impliedCopyNumber;
        }

        double purity = purityAdjuster.purity();
        double normFactor = purityAdjuster.normFactor();

        // Minimise
        final double minBAFTotalDeviation =
                mPloidyDeviation.majorAlleleDeviation(purity, normFactor, minBAFMajorAllelePloidy) + mPloidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        minBAFMinorAllelePloidy);
        final double maxBAFTotalDeviation =
                mPloidyDeviation.majorAlleleDeviation(purity, normFactor, maxBAFMajorAllelePloidy) + mPloidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        maxBAFMinorAllelePloidy);
        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? 0.5 : mAmbiguousBaf;
    }
}
