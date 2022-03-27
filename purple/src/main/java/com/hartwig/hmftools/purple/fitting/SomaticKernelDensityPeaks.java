package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.purple.somatic.SomaticData;

import org.jetbrains.annotations.NotNull;

final class SomaticKernelDensityPeaks
{
    private static final double KERNEL_BANDWIDTH = 0.03;

    public static Optional<FittedPurity> fitPurity(
            final List<FittedPurity> allCandidates, final List<SomaticData> variants,
            int minVariants, double minPeak, double minPurity, double maxPurity)
    {
        if(variants.size() < minVariants)
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", variants.size());
            return Optional.empty();
        }

        final List<SomaticPeak> peaks = SomaticKernelDensityPeaks.findSomaticPeaks(variants);

        // First try and get the largest implied purity where peak count > minPeak
        int maxPeak = 0;
        for(int i = peaks.size() - 1; i >= 0; i--)
        {
            SomaticPeak peak = peaks.get(i);
            double impliedPurity = peak.alleleFrequency() * 2;
            if(inPurityRange(impliedPurity, minPurity, maxPurity))
            {
                if(peak.count() >= minPeak)
                {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if(diploid.isPresent())
                    {
                        PPL_LOGGER.debug("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    }
                    else
                    {
                        PPL_LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
                maxPeak = Math.max(maxPeak, peak.count());
            }
        }

        // Failing that, get the implied purity with the largest peak
        if(maxPeak > 0)
        {
            for(int i = peaks.size() - 1; i >= 0; i--)
            {
                SomaticPeak peak = peaks.get(i);
                double impliedPurity = peak.alleleFrequency() * 2;
                if(peak.count() == maxPeak)
                {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if(diploid.isPresent())
                    {
                        PPL_LOGGER.debug("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    }
                    else
                    {
                        PPL_LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
            }
        }

        PPL_LOGGER.debug("Unable to determine somatic implied purity.");
        return Optional.empty();
    }

    private static boolean inPurityRange(double impliedPurity, double minPurity, double maxPurity)
    {
        return Doubles.greaterOrEqual(impliedPurity, minPurity) && Doubles.lessOrEqual(impliedPurity, maxPurity);
    }

    private static Optional<FittedPurity> diploid(double purity, List<FittedPurity> diploidCandidates)
    {
        return diploidCandidates.stream().filter(x -> Doubles.equal(x.purity(), purity)).findFirst();
    }

    public static List<SomaticPeak> findSomaticPeaks(final List<SomaticData> variants)
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
                final SomaticPeak peak = ImmutableSomaticPeak.builder().alleleFrequency(alleleFrequency).count(peakCount).build();
                PPL_LOGGER.debug("discovered peak {}", peak);
                results.add(peak);
            }
        }

        return results;
    }

    private static int count(double peak, @NotNull final List<Double> sample)
    {
        return (int) sample.stream().filter(vaf -> between(peak, vaf - 0.015, vaf + 0.015)).count();
    }

    private static boolean between(double victim, double min, double max)
    {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }
}
