package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.jetbrains.annotations.Nullable;

public final class SomaticKernelDensityPeaks
{
    private static final double KERNEL_BANDWIDTH = 0.03;

    public static @Nullable FittedPurity fitPurity(
            final List<FittedPurity> allCandidates, final List<SomaticVariant> variants,
            int minVariants, double minPeak, double minPurity, double maxPurity)
    {
        if(variants.size() < minVariants)
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", variants.size());
            return null;
        }

        final List<SomaticPeak> peaks = SomaticKernelDensityPeaks.findSomaticPeaks(variants);

        // First try and get the largest implied purity where peak count > minPeak
        int maxPeak = 0;
        for(int i = peaks.size() - 1; i >= 0; i--)
        {
            SomaticPeak peak = peaks.get(i);
            double impliedPurity = peak.AlleleFrequency * 2;
            if(inPurityRange(impliedPurity, minPurity, maxPurity))
            {
                if(peak.Count >= minPeak)
                {
                    FittedPurity matchedFittedPurity = findMatchedFittedPurity(impliedPurity, allCandidates, 0.001);
                    if(matchedFittedPurity != null)
                    {
                        PPL_LOGGER.debug("somatic implied purity({})", impliedPurity);
                        return matchedFittedPurity;
                    }
                    else
                    {
                        PPL_LOGGER.warn("unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }

                maxPeak = Math.max(maxPeak, peak.Count);
            }
        }

        // Failing that, get the implied purity with the largest peak
        if(maxPeak > 0)
        {
            for(int i = peaks.size() - 1; i >= 0; i--)
            {
                SomaticPeak peak = peaks.get(i);
                double impliedPurity = peak.AlleleFrequency * 2;
                if(peak.Count == maxPeak)
                {
                    FittedPurity matchedFittedPurity = findMatchedFittedPurity(impliedPurity, allCandidates, 0.001);
                    if(matchedFittedPurity != null)
                    {
                        PPL_LOGGER.debug("somatic implied purity({})", impliedPurity);
                        return matchedFittedPurity;
                    }
                    else
                    {
                        PPL_LOGGER.warn("unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
            }
        }

        PPL_LOGGER.debug("unable to determine somatic implied purity.");
        return null;
    }

    private static boolean inPurityRange(double impliedPurity, double minPurity, double maxPurity)
    {
        return Doubles.greaterOrEqual(impliedPurity, minPurity) && Doubles.lessOrEqual(impliedPurity, maxPurity);
    }

    protected static FittedPurity findMatchedFittedPurity(double purity, final List<FittedPurity> allCandidates, final double epsilon)
    {
        return allCandidates.stream().filter(x -> abs(x.purity() - purity) < epsilon).findFirst().orElse(null);
    }

    public static List<SomaticPeak> findSomaticPeaks(final List<SomaticVariant> variants)
    {
        return findPeaks(variants.stream().map(x -> x.alleleFrequency()).collect(Collectors.toList()));
    }

    public static List<SomaticPeak> findPeaks(final List<Double> sample)
    {
        final KernelEstimator estimator = new KernelEstimator(0.001, KERNEL_BANDWIDTH);
        sample.forEach(x -> estimator.addValue(x, 1.0D));

        final double[] vafs = IntStream.rangeClosed(0, 51).mapToDouble(x -> x / 100d).toArray();
        final double[] densities = DoubleStream.of(vafs).map(estimator::getProbability).toArray();

        final List<SomaticPeak> results = Lists.newArrayList();
        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];
            if(Doubles.greaterThan(density, densities[i - 1]) && Doubles.greaterThan(density, densities[i + 1]))
            {
                final double alleleFrequency = vafs[i];
                final int peakCount = count(alleleFrequency, sample);

                PPL_LOGGER.debug(format("discovered somatic peak: count(%d) alleleFrequency(%.3f)", peakCount, alleleFrequency));
                results.add(new SomaticPeak(alleleFrequency, peakCount));
            }
        }

        return results;
    }

    private static int count(double peak, final List<Double> sample)
    {
        return (int) sample.stream().filter(vaf -> between(peak, vaf - 0.015, vaf + 0.015)).count();
    }

    private static boolean between(double victim, double min, double max)
    {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }
}
