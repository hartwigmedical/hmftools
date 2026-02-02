package com.hartwig.hmftools.amber.contamination;

import static java.lang.Math.floor;
import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.MIN_NORMAL_READ_DEPTH;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_MIN;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_SITE_PERC;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_SITE_LOW_VAF_PERC;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_VAF_MIN;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.Integers;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class TumorContaminationModel
{
    private static final double INCREMENT = 0.001;

    public double calcContamination(final List<TumorContamination> unfilteredTumorContamination, long amberSiteCount)
    {
        //        List<TumorContamination> contaminationSites = unfilteredTumorContamination.stream()
        //                .filter(x -> x.Normal.readDepth() > MIN_NORMAL_READ_DEPTH).collect(Collectors.toList());

        List<TumorContamination> contaminationSites = unfilteredTumorContamination.stream()
                .filter(x -> x.Tumor.readDepth() > MIN_NORMAL_READ_DEPTH)
                .filter(x ->
                {
                    double v = x.tumorVaf();
                    return v > 0.001 && v < 0.2;
                })
                .collect(Collectors.toList());

        int medianTumorReadDepth = medianDepth(contaminationSites);

        AMB_LOGGER.debug("contamination sites median tumor depth({})", medianTumorReadDepth);

        long threePlusReadsSiteCount = 0;
        long threePlusReadsLowVafSiteCount = 0;
        long twoPlusReadsSiteCount = 0;

        for(TumorContamination site : contaminationSites)
        {
            if(site.Tumor.altSupport() >= 2)
            {
                ++twoPlusReadsSiteCount;

                if(site.Tumor.altSupport() >= 3)
                {
                    ++threePlusReadsSiteCount;

                    double vaf = site.Tumor.altSupport() / (double) site.Tumor.readDepth();
                    if(vaf < THREE_PLUS_READS_VAF_MIN)
                    {
                        ++threePlusReadsLowVafSiteCount;
                    }
                }
            }
        }

        boolean calcContamination = false;

        if(threePlusReadsSiteCount >= THREE_PLUS_READS_MIN)
        {
            if(threePlusReadsSiteCount > THREE_PLUS_READS_SITE_PERC * amberSiteCount)
            {
                calcContamination = true;
            }
            else if(threePlusReadsLowVafSiteCount > THREE_PLUS_READS_SITE_LOW_VAF_PERC * amberSiteCount)
            {
                calcContamination = true;
            }
        }

        if(!calcContamination)
        {
            AMB_LOGGER.debug("siteCount({}) altSites(3={} 2={} 3-lowVaf={})",
                    amberSiteCount, threePlusReadsSiteCount, twoPlusReadsSiteCount, threePlusReadsLowVafSiteCount);
            return 0;
        }

        Map<Integer, Long> altSupportFrequencies = contaminationSites.stream()
                .collect(Collectors.groupingBy(x -> x.Tumor.altSupport(), Collectors.counting()));

        double contaminationLevel = doCalcContamination(medianTumorReadDepth, twoPlusReadsSiteCount, altSupportFrequencies);
        AMB_LOGGER.info("contamination level identified({}) ", floor(contaminationLevel * 1000) / 10);

        AMB_LOGGER.debug("contamination({}) siteCount({}) altSites(3={} 2={} 3-lowVaf={})",
                format("%.6f", contaminationLevel), amberSiteCount, threePlusReadsSiteCount, twoPlusReadsSiteCount, threePlusReadsLowVafSiteCount);

        return contaminationLevel;
    }

    @VisibleForTesting
    public double calcContamination(int medianTumorReadDepth, final Map<Integer, Long> altSupportFrequencies)
    {
        long twoPlusReadCount = calcAltReadCountAboveThreshold(2, altSupportFrequencies);
        return doCalcContamination(medianTumorReadDepth, twoPlusReadCount, altSupportFrequencies);
    }

    private double doCalcContamination(
            int medianTumorReadDepth, long twoPlusReadCount, final Map<Integer, Long> altSupportFrequencies)
    {
        double contamination = 0;
        double lowestScore = Double.MAX_VALUE;

        for(double i = INCREMENT; Doubles.lessOrEqual(i, 1); i = i + INCREMENT)
        {
            double score = contaminationScore(i, twoPlusReadCount, medianTumorReadDepth, altSupportFrequencies);
            if(Doubles.lessThan(score, lowestScore))
            {
                lowestScore = score;
                contamination = i;
            }
        }

        return contamination;
    }

    private double contaminationScore(
            final double contamination, final long twoPlusReadCount, final long medianTumorReadDepth,
            final Map<Integer, Long> altSupportMap)
    {
        PoissonDistribution hetDistribution = new PoissonDistribution(0.5 * contamination * medianTumorReadDepth);
        PoissonDistribution homAltDistribution = new PoissonDistribution(contamination * medianTumorReadDepth);

        Function<Integer, Double> unadjustedModelLikelihood = altSupport ->
        {
            double hetLikelihood = hetDistribution.probability(altSupport);
            double homAltLikelihood = homAltDistribution.probability(altSupport);
            return 0.5 * hetLikelihood + 0.25 * homAltLikelihood;
        };

        double altSupport0 = unadjustedModelLikelihood.apply(0) + 0.25 * 1;
        double altSupport1 = unadjustedModelLikelihood.apply(1);
        double altSupport2PlusAdjustment = 1 - altSupport0 - altSupport1;

        Function<Integer, Double> modelLikelihood =
                integer -> unadjustedModelLikelihood.apply(integer) / altSupport2PlusAdjustment;

        double totalModelPercentage = 0;
        double totalDifference = 0;

        for(Map.Entry<Integer, Long> entry : altSupportMap.entrySet())
        {
            long altSupport = entry.getKey();
            long altSupportCount = entry.getValue();

            if(altSupport > 1)
            {
                double modelPercentage = modelLikelihood.apply((int) altSupport);
                double actualPercentage = altSupportCount * 1d / twoPlusReadCount;

                totalDifference += Math.abs(actualPercentage - modelPercentage);
                totalModelPercentage += modelPercentage;
            }
        }

        return totalDifference + (1 - totalModelPercentage);
    }

    static long calcAltReadCountAboveThreshold(int minAltSupport, final Map<Integer, Long> altSupportMap)
    {
        return altSupportMap.entrySet().stream().filter(x -> x.getKey() >= minAltSupport).mapToLong(Map.Entry::getValue).sum();
    }

    @VisibleForTesting
    static int medianDepth(final List<TumorContamination> contamination)
    {
        return Integers.medianPositiveValue(contamination.stream().map(x -> x.Tumor.readDepth()).collect(Collectors.toList()));
    }
}
