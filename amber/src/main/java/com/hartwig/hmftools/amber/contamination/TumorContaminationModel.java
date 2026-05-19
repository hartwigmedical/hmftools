package com.hartwig.hmftools.amber.contamination;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.MIN_NORMAL_READ_DEPTH;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_MIN;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_SITE_PERC;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_SITE_LOW_VAF_PERC;
import static com.hartwig.hmftools.amber.AmberConstants.THREE_PLUS_READS_VAF_MIN;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.Integers;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class TumorContaminationModel
{
    private static final double INCREMENT = 0.001;

    public double calcContamination(final List<TumorContamination> allTumorContamination, long amberSiteCount)
    {
        List<TumorContamination> contaminationSites = allTumorContamination.stream()
                .filter(x -> x.Normal.readDepth() > MIN_NORMAL_READ_DEPTH).collect(Collectors.toList());

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

        int medianTumorReadDepth =  Integers.medianPositiveValue(contaminationSites.stream()
                .map(x -> x.Tumor.readDepth()).collect(Collectors.toList()));

        double contaminationLevel = calcContamination(medianTumorReadDepth, twoPlusReadsSiteCount, altSupportFrequencies);

        AMB_LOGGER.info("contamination level identified({}) ", format("%.3f", contaminationLevel));

        AMB_LOGGER.debug("contamination({}) siteCount({}) medianTumorDepth({}) altSites(3={} 2={} 3-lowVaf={})",
                format("%.3f", contaminationLevel), amberSiteCount, medianTumorReadDepth,
                threePlusReadsSiteCount, twoPlusReadsSiteCount, threePlusReadsLowVafSiteCount);

        return contaminationLevel;
    }

    @VisibleForTesting
    public double calcContamination(int medianTumorReadDepth, long twoPlusReadCount, final Map<Integer,Long> altSupportFrequencies)
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

        double altSupport0 = calcUnadjustedLikelihood(hetDistribution, homAltDistribution, 0) + 0.25;
        double altSupport1 = calcUnadjustedLikelihood(hetDistribution, homAltDistribution, 1);
        double altSupport2Plus = 1 - altSupport0 - altSupport1;

        double totalModelPercentage = 0;
        double totalDifference = 0;

        for(Map.Entry<Integer, Long> entry : altSupportMap.entrySet())
        {
            long altSupport = entry.getKey();
            long altSupportCount = entry.getValue();

            if(altSupport > 1)
            {
                double modelPercentage = calcUnadjustedLikelihood(hetDistribution, homAltDistribution, (int)altSupport) / altSupport2Plus;

                double actualPercentage = altSupportCount * 1d / twoPlusReadCount;

                totalDifference += Math.abs(actualPercentage - modelPercentage);
                totalModelPercentage += modelPercentage;
            }
        }

        return totalDifference + (1 - totalModelPercentage);
    }

    private static double calcUnadjustedLikelihood(
            final PoissonDistribution hetDistribution, final PoissonDistribution homAltDistribution, int altSupport)
    {
        double hetLikelihood = hetDistribution.probability(altSupport);
        double homAltLikelihood = homAltDistribution.probability(altSupport);
        return 0.5 * hetLikelihood + 0.25 * homAltLikelihood;
    }
}
