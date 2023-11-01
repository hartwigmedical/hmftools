package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityResult.INVALID_RESULT;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.wisp.common.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

public class VafPeakModel extends ClonalityModel
{
    private double mTumorAvgVaf;

    public VafPeakModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        super(config, resultsWriter, sample,  variants);

        mTumorAvgVaf = calcTumorAvgVaf();
    }

    @Override
    public ClonalityResult calculate(final String sampleId, final FragmentCalcResult estimatedResult)
    {
        if(estimatedResult.PurityProbability > PurityConstants.SOMATIC_PEAK_MAX_PROBABILITY)
            return INVALID_RESULT;

        List<SomaticVariant> filteredVariants = Lists.newArrayList();
        List<VariantCalcs> variantCalcData = Lists.newArrayList();

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null)
                continue;

            if(canUseVariant(variant, sampleFragData))
            {
                filteredVariants.add(variant);

                VariantCalcs varCalcData = calcAdjustedVaf(sampleFragData, tumorFragData);
                variantCalcData.add(varCalcData);
            }
        }

        if(filteredVariants.size() < PurityConstants.SOMATIC_PEAK_MIN_VARIANTS)
            return INVALID_RESULT;

        List<Double> variantVafs = variantCalcData.stream().map(x -> x.SampleVaf).collect(Collectors.toList());

        List<VafPeak> vafPeaks = findVafPeaks(variantVafs, estimatedResult.VAF);

        if(!vafPeaks.isEmpty())
        {
            VafPeak maxVafPeak = vafPeaks.get(vafPeaks.size() - 1);
            VafPeak minVafPeak = vafPeaks.get(0);

            CT_LOGGER.debug("sample({}) filteredVars({}) found {} somatic vaf peaks, max({})",
                    sampleId, filteredVariants.size(), vafPeaks.size(), format("%.3f", maxVafPeak.Peak));

            for(VafPeak vafPeak : vafPeaks)
            {
                CT_LOGGER.debug("sample({}) somatic vaf peak({})", sampleId, vafPeak);
            }

            return new ClonalityResult(ClonalityMethod.VAF_PEAK, maxVafPeak.Peak, maxVafPeak.Peak, minVafPeak.Peak, maxVafPeak.Count, 0);
        }

        return INVALID_RESULT;
    }

    private class VafPeak implements Comparable<VafPeak>
    {
        public final double Peak;
        public final int Count;

        public VafPeak(final double peak, final int count)
        {
            Peak = peak;
            Count = count;
        }

        @Override
        public int compareTo(final VafPeak other)
        {
            if(Peak == other.Peak)
                return 0;

            return Peak < other.Peak ? -1 : 1;
        }

        public String toString() { return format("%.3f=%d", Peak, Count); }
    }

    private static final double PEAK_VAF_BUFFER = 0.015;

    private List<VafPeak> findVafPeaks(final List<Double> sampleVafs, double estimatedVaf)
    {
        double vafTotal = 0;
        double maxVaf = 0;
        for(double variantVaf : sampleVafs)
        {
            maxVaf = max(maxVaf, variantVaf);
            vafTotal += variantVaf;
        }

        if(vafTotal == 0)
            return Collections.emptyList();

        double avgVaf = vafTotal / sampleVafs.size();

        // double densityBandwidth = min(avgVaf / 2, 0.01);
        double densityBandwidth = max(avgVaf/8, min(avgVaf/2, 0.01));
        // bandwidth=pmax(mean(filteredVars$SampleVaf)/8,pmin(mean(filteredVars$SampleVaf)/2,0.01))  # mean = 0.076, min(avgVaf / 2, 0.01)
        int maxVafLimit = min((int)round(maxVaf * 100), 99);

        int vafFraction = 5;
        double[] vafs = IntStream.rangeClosed(0, maxVafLimit * vafFraction).mapToDouble(x -> x / (100d * vafFraction)).toArray();

        KernelEstimator estimator = new KernelEstimator(0.001, densityBandwidth);

        sampleVafs.forEach(x -> estimator.addValue(x, 1.0));

        double[] densities = DoubleStream.of(vafs).map(estimator::getProbability).toArray();

        final List<VafPeak> peakVafs = Lists.newArrayList();

        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];

            // peak must be above the estimated VAF

            if(!Doubles.greaterThan(density, densities[i - 1]) || !Doubles.greaterThan(density, densities[i + 1]))
                continue;

            double densityVaf = vafs[i];

            if(densityVaf < estimatedVaf)
                continue;

            // count up observations at this density peak
            int peakCount = 0;

            for(Double variantVar : sampleVafs)
            {
                if(variantVar >= densityVaf - PEAK_VAF_BUFFER && variantVar <= densityVaf + PEAK_VAF_BUFFER)
                    ++peakCount;
            }

            if(peakCount < PurityConstants.SOMATIC_PEAK_MIN_PEAK_VARIANTS)
                continue;

            CT_LOGGER.debug(format("somatic peak: count(%d) vaf(%.3f) densityBandwidth(%.4f)", peakCount, densityVaf, densityBandwidth));
            peakVafs.add(new VafPeak(densityVaf, peakCount));
        }

        Collections.sort(peakVafs);

        return peakVafs;
    }

    private boolean canUseVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return useVariant(variant, sampleFragData)
                && sampleFragData.UmiCounts.totalCount() >= PurityConstants.SOMATIC_PEAK_MIN_DEPTH
                && sampleFragData.UmiCounts.alleleCount() >= PurityConstants.SOMATIC_PEAK_MIN_AD;
    }

    private class VariantCalcs
    {
        public final double TumorVaf;
        public final double TumorVafAdjust;
        public final double SampleVaf;

        public VariantCalcs(final double tumorVaf, final double tumorVafAdjust, final double sampleVaf)
        {
            TumorVaf = tumorVaf;
            TumorVafAdjust = tumorVafAdjust;
            SampleVaf = sampleVaf;
        }

        public double sampleAdjustedVaf() { return SampleVaf * TumorVafAdjust; }
    }

    private VariantCalcs calcAdjustedVaf(final GenotypeFragments sampleFragData, final GenotypeFragments tumorFragData)
    {
        double variantVaf = sampleFragData.vaf();
        double tumorVaf = tumorFragData.vaf();

        if(tumorVaf == 0)
            return new VariantCalcs(0, 0, 0);

        double tumorVafAdjust = mTumorAvgVaf / tumorVaf;
        return new VariantCalcs(tumorVaf, tumorVafAdjust, variantVaf);
    }

    private double calcTumorAvgVaf()
    {
        double sampleVafTotal = 0;
        int sampleVafCount = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(tumorFragData == null)
                continue;

            if(!variant.isFiltered())
            {
                if(tumorFragData.Depth > 0)
                {
                    sampleVafTotal += tumorFragData.vaf();
                    ++sampleVafCount;
                }
            }
        }

        return sampleVafCount > 0 ? sampleVafTotal / sampleVafCount : 0;
    }
}
