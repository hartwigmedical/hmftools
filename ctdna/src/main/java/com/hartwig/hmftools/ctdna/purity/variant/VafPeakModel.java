package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_MIN_DEPTH;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_PROBABILITY;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_VAF_INCREMENT;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SOMATIC_PEAK_MAX_PROBABILITY;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SOMATIC_PEAK_MIN_AD;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SOMATIC_PEAK_MIN_DEPTH;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SOMATIC_PEAK_MIN_PEAK_VARIANTS;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SOMATIC_PEAK_MIN_VARIANTS;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.DROPOUT_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.ctdna.common.SampleData;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.ResultsWriter;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.checkerframework.checker.units.qual.C;

public class VafPeakModel
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<SomaticVariant> mVariants;

    // calculate state
    private double mTumorAvgVaf;
    private double mMaxSomaticVaf;

    public VafPeakModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mVariants = variants;

        mTumorAvgVaf = calcTumorAvgVaf();
        mMaxSomaticVaf = 0;
    }

    public double maxSomaticVaf() { return mMaxSomaticVaf; }

    public void calculate(final String sampleId, final FragmentCalcResult estimatedResult)
    {
        if(estimatedResult.PurityProbability > SOMATIC_PEAK_MAX_PROBABILITY)
            return;

        List<SomaticVariant> filteredVariants = Lists.newArrayList();
        List<VariantCalcs> variantCalcData = Lists.newArrayList();

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null)
                continue;

            if(useVariant(variant, sampleFragData))
            {
                filteredVariants.add(variant);

                VariantCalcs varCalcData = calcAdjustedVaf(sampleFragData, tumorFragData);
                variantCalcData.add(varCalcData);

                if(varCalcData.SampleVaf > 0)
                {
                    double sampleAdjVaf = varCalcData.sampleAdjustedVaf();
                    // maxSampleVaf = max(maxSampleVaf, sampleAdjVaf);
                    // minSampleVaf = minSampleVaf == 0 ? sampleAdjVaf : min(minSampleVaf, sampleAdjVaf);
                }
            }
        }

        if(filteredVariants.size() < SOMATIC_PEAK_MIN_VARIANTS)
            return;

        List<Double> variantVafs = variantCalcData.stream().map(x -> x.SampleVaf).collect(Collectors.toList());

        List<VafPeak> vafPeaks = findVafPeaks(variantVafs, estimatedResult.VAF);

        if(!vafPeaks.isEmpty())
        {
            mMaxSomaticVaf = vafPeaks.stream().mapToDouble(x -> x.Peak).max().orElse(0);

            CT_LOGGER.debug("sample({}) found {} somatic vaf peaks, max({})",
                    sampleId, vafPeaks.size(), format("%.3f", mMaxSomaticVaf));
        }
    }

    private class VafPeak
    {
        public final double Peak;
        public final int Count;

        public VafPeak(final double peak, final int count)
        {
            Peak = peak;
            Count = count;
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

        double densityBandwidth = min(avgVaf / 2, 0.01);
        int maxVafLimit = min((int)round(maxVaf * 100), 99);

        int vafFraction = 5;
        double[] vafs = IntStream.rangeClosed(0, maxVafLimit * vafFraction).mapToDouble(x -> x / (100d * vafFraction)).toArray();

        KernelEstimator estimator = new KernelEstimator(0.001, vafs.length * 2);

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

            if(peakCount < SOMATIC_PEAK_MIN_PEAK_VARIANTS)
                continue;

            CT_LOGGER.debug(format("somatic peak: count(%d) vaf(%.3f)", peakCount, densityVaf));
            peakVafs.add(new VafPeak(densityVaf, peakCount));
        }

        return peakVafs;
    }

    private boolean useVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return variant.PassFilters
                && variant.sequenceGcRatio() >= mConfig.GcRatioMin
                && sampleFragData.UmiCounts.total() >= SOMATIC_PEAK_MIN_DEPTH
                && sampleFragData.UmiCounts.alleleTotal() >= SOMATIC_PEAK_MIN_AD;
    }

    public class VariantSimVafCalcs
    {
        public final double ExpectedAD;
        public final double Probability;

        public VariantSimVafCalcs(final double expectedAD, final double probability)
        {
            ExpectedAD = expectedAD;
            Probability = probability;
        }
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

            if(variant.PassFilters)
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

    public static BufferedWriter initialiseWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(DROPOUT_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(config.multipleSamples())
                writer.write("PatientId\t");

            writer.write("SampleId\tSimulatedVaf\tVariantCount\tDropoutCount");
            writer.write("\tProbabilityScore\tLowMeanDiffTotal\tHighMeanDiffTotal\tVafStdDev\tZScoreAvg");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    /*
    public synchronized static void writeCalcData(
            final BufferedWriter writer, final PurityConfig config, final SampleData sample, final String sampleId,
            int varCount, final SimulatedVafCalcs simVafCalcs)
    {
        if(writer == null)
            return;

        try
        {
            if(config.multipleSamples())
                writer.write(format("%s\t", sample.PatientId));

            writer.write(format("%s\t%.2f\t%d\t%d", sampleId, simVafCalcs.SimulatedVaf, varCount, simVafCalcs.DropoutCount));

            writer.write(format("\t%4.3f\t%.3f\t%.3f\t%.3f\t%.3f",
                    simVafCalcs.ProbabilityScore, simVafCalcs.LowMeanDiffTotal, simVafCalcs.HighMeanDiffTotal,
                    simVafCalcs.VafStdDev, simVafCalcs.ZScoreAvg));
            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write dropout calc file: {}", e.toString());
        }

    }
    */
}
