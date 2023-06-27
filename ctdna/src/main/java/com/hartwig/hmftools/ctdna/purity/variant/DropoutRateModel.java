package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_MIN_DEPTH;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_PROBABILITY;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_VAF_INCREMENT;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.DROPOUT_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.ResultsWriter;
import com.hartwig.hmftools.ctdna.purity.SampleData;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class DropoutRateModel
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<SomaticVariant> mVariants;

    // calculate state
    private double mTumorAvgVaf;

    public DropoutRateModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mVariants = variants;

        mTumorAvgVaf = calcTumorAvgVaf();
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

    private class SimulatedVafCalcs
    {
        public final double SimulatedVaf;
        public final int DropoutCount;
        public final double ProbabilityScore;
        public final double LowMeanDiffTotal;
        public final double HighMeanDiffTotal;
        public final double VafStdDev;
        public final double ZScoreAvg;

        public SimulatedVafCalcs(
                final double simulatedVaf, final int dropoutCount, final double probabilityScore,
                final double lowMeanDiffTotal, final double highMeanDiffTotal, final double vafStdDev, final double zScoreAvg)
        {
            SimulatedVaf = simulatedVaf;
            DropoutCount = dropoutCount;
            ProbabilityScore = probabilityScore;
            ZScoreAvg = zScoreAvg;
            LowMeanDiffTotal = lowMeanDiffTotal;
            HighMeanDiffTotal = highMeanDiffTotal;
            VafStdDev = vafStdDev;
        }
    }

    private static final double MAX_LOG_PROB = abs(log10(DROPOUT_RATE_PROBABILITY));
    private static final double MAX_Z_SCORE = 4;

    public void calculate(final String sampleId)
    {
        List<SomaticVariant> filteredVariants = Lists.newArrayList();
        List<VariantCalcs> variantCalcData = Lists.newArrayList();

        double minSampleVaf = 0;
        double maxSampleVaf = 0;

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
                    maxSampleVaf = max(maxSampleVaf, sampleAdjVaf);
                    minSampleVaf = minSampleVaf == 0 ? sampleAdjVaf : min(minSampleVaf, sampleAdjVaf);
                }
            }
        }

        double simVafIncrement = DROPOUT_RATE_VAF_INCREMENT; // may set from config
        double minSimulatedVaf = ceil(minSampleVaf / simVafIncrement) * simVafIncrement;
        double maxSimulatedVaf = floor(maxSampleVaf / simVafIncrement) * simVafIncrement;
        int varCount = variantCalcData.size();

        List<SimulatedVafCalcs> simulatedVafCalcs = Lists.newArrayList();

        for(double simulatedVaf = minSimulatedVaf; simulatedVaf <= maxSimulatedVaf; simulatedVaf += simVafIncrement)
        {
            List<VariantSimVafCalcs> simVarCalcData = Lists.newArrayListWithCapacity(varCount);

            double probabilityScore = 0;
            double meanDiffTotal = 0;
            int dropoutCount = 0;

            List<Double> meanVafDiffs = Lists.newArrayList();

            // calculate a standard deviation around this simulated mean
            for(int i = 0; i < varCount; ++i)
            {
                SomaticVariant variant = filteredVariants.get(i);
                GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
                VariantCalcs varCalcData = variantCalcData.get(i);

                VariantSimVafCalcs varSimCalcData = calcSimulatedSampleData(simulatedVaf, sampleFragData, varCalcData);

                boolean belowExpected = sampleFragData.AlleleCount < varSimCalcData.ExpectedAD;

                if(belowExpected && varSimCalcData.Probability < DROPOUT_RATE_PROBABILITY)
                {
                    ++dropoutCount;
                    continue;
                }

                double logProb = varSimCalcData.Probability > 0 ? log10(varSimCalcData.Probability) : MAX_LOG_PROB;

                probabilityScore += max(MAX_LOG_PROB - abs(logProb), 0);

                double vafDiff = varCalcData.sampleAdjustedVaf() - simulatedVaf;
                meanVafDiffs.add(vafDiff);
                meanDiffTotal += pow(vafDiff, 2);
            }

            double simVafStdDev = sqrt(meanDiffTotal / (varCount - 1));
            double lowMeanDiffTotal = 0;
            double highMeanDiffTotal = 0;
            double zScoreTotal = 0;

            for(Double vafDiff : meanVafDiffs)
            {
                if(vafDiff > 0)
                    highMeanDiffTotal += vafDiff;
                else
                    lowMeanDiffTotal += abs(vafDiff);

                zScoreTotal += max(MAX_Z_SCORE - abs(vafDiff), 0);
            }

            double zScoreAvg = zScoreTotal > 0 ? zScoreTotal / (varCount - dropoutCount) : 0;

            SimulatedVafCalcs simVafCalcData = new SimulatedVafCalcs(
                simulatedVaf, dropoutCount, probabilityScore, lowMeanDiffTotal, highMeanDiffTotal, simVafStdDev, zScoreAvg);

            simulatedVafCalcs.add(simVafCalcData);

            writeCalcData(mResultsWriter.getDropoutWriter(), mConfig, mSample, sampleId, varCount, simVafCalcData);

            // compute cost functions for the simulated results

        }
    }

    private boolean useVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return variant.PassFilters
                && variant.sequenceGcRatio() >= mConfig.GcRatioMin
                && sampleFragData.UmiCounts.total() >= DROPOUT_RATE_MIN_DEPTH;
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

    private VariantSimVafCalcs calcSimulatedSampleData(
            double simulatedVaf, final GenotypeFragments sampleFragData, final VariantCalcs variantCalcs)
    {
        double expectedAD = simulatedVaf * variantCalcs.TumorVafAdjust * sampleFragData.Depth;
        int actualAD = sampleFragData.UmiCounts.alleleTotal();

        double probability;
        PoissonDistribution poisson = new PoissonDistribution(expectedAD);

        if(actualAD > expectedAD)
            probability = 1 - poisson.cumulativeProbability(actualAD - 1);
        else
            probability = poisson.cumulativeProbability(actualAD);

        return new VariantSimVafCalcs(expectedAD, probability);
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

}
