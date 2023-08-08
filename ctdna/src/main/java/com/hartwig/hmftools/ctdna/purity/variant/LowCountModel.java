package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DROPOUT_RATE_INCREMENT;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.DROPOUT_FILE_ID;
import static com.hartwig.hmftools.ctdna.purity.variant.ClonalityResult.INVALID_RESULT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.ResultsWriter;
import com.hartwig.hmftools.ctdna.common.SampleData;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class LowCountModel extends ClonalityModel
{
    // calculate state
    private double mTumorAvgVaf;

    public LowCountModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        super(config, resultsWriter, sample,  variants);

        mTumorAvgVaf = calcTumorAvgVaf();
    }

    @Override
    public ClonalityResult calculate(final String sampleId, final FragmentCalcResult estimatedResult)
    {
        if(estimatedResult.VAF == 0)
            return INVALID_RESULT;

        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        int observedFrag1 = 0;
        int observedFrag2Plus = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(useVariant(variant, sampleFragData) && sampleFragData.Depth > 0)
            {
                filteredVariants.add(variant);

                if(sampleFragData.AlleleCount == 1)
                    ++observedFrag1;
                else if(sampleFragData.AlleleCount == 2)
                    ++observedFrag2Plus;
            }
        }

        if(observedFrag1 == 0)
            return INVALID_RESULT;

        List<SimulatedVafCalcs> simulatedVafCalcs = Lists.newArrayList();

        // now test each simluated dropout rate and VAF
        for(double dropoutRate = DROPOUT_RATE_INCREMENT; dropoutRate < 0.95; dropoutRate += DROPOUT_RATE_INCREMENT)
        {
            double simulatedVaf = estimatedResult.VAF / (1 - dropoutRate);

            if(simulatedVaf >= 1)
                break;

            double probTotalFrag1 = 0;
            double probTotalFrag2Plus = 0;

            for(SomaticVariant variant : filteredVariants)
            {
                GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

                double expectedAD = sampleFragData.Depth * simulatedVaf;

                PoissonDistribution poisson = new PoissonDistribution(expectedAD);

                double probFrag1 = poisson.probability(1);
                double probFrag2Plus = 1 - poisson.cumulativeProbability(1);

                probTotalFrag1 += probFrag1;
                probTotalFrag2Plus += probFrag2Plus;
            }

            SimulatedVafCalcs simVafCalcs = new SimulatedVafCalcs(dropoutRate, simulatedVaf, probTotalFrag1, probTotalFrag2Plus);
            simulatedVafCalcs.add(simVafCalcs);

            writeSimulatedDropoutData(
                    mResultsWriter.getDropoutWriter(), mConfig, mSample, sampleId, filteredVariants.size(), simVafCalcs,
                    observedFrag1, observedFrag2Plus);
        }

        // now find the closest ratio to the observed ratio
        double observedRatio = observedFrag2Plus / (double)observedFrag1;
        double calcVaf = findVafRatio(observedRatio, simulatedVafCalcs);

        double lowObservedRatio = max(observedFrag2Plus - 2, 1) / (double)observedFrag1;
        double lowCalcVaf = findVafRatio(lowObservedRatio, simulatedVafCalcs);

        double highObservedRatio = (observedFrag2Plus + 2) / (double)observedFrag1;
        double highCalcVaf = findVafRatio(highObservedRatio, simulatedVafCalcs);

        CT_LOGGER.debug(format("sample(%s) low-count model: obsRatio(%.2f 1=%d 2+=%d) vaf((%.6f low=%.6f high=%.6f)",
                sampleId, observedRatio, observedFrag1, observedFrag2Plus, calcVaf, lowCalcVaf, highCalcVaf));

        return new ClonalityResult(ClonalityMethod.LOW_COUNT, calcVaf, lowCalcVaf, highCalcVaf, observedFrag1 + observedFrag2Plus);
    }

    private static double findVafRatio(double observedRatio, final List<SimulatedVafCalcs> simulatedVafCalcs)
    {
        SimulatedVafCalcs closestRatioUp = null;
        SimulatedVafCalcs closestRatioDown = null;

        for(SimulatedVafCalcs simVafCalcs : simulatedVafCalcs)
        {
            double probRatio = simVafCalcs.fragmentRatio();

            if(probRatio > observedRatio)
            {
                if(closestRatioUp == null || probRatio < closestRatioUp.fragmentRatio())
                    closestRatioUp = simVafCalcs;
            }
            else
            {
                if(closestRatioDown == null || probRatio > closestRatioDown.fragmentRatio())
                    closestRatioDown = simVafCalcs;
            }
        }

        double calcVaf = 0;

        if(closestRatioUp != null && closestRatioDown != null)
        {
            double upperFraction = (observedRatio - closestRatioDown.fragmentRatio())
                    / (closestRatioUp.fragmentRatio() - closestRatioDown.fragmentRatio());

            calcVaf = upperFraction * closestRatioUp.SimulatedVaf + (1 - upperFraction) * closestRatioDown.SimulatedVaf;
        }
        else if(closestRatioUp != null)
        {
            calcVaf = closestRatioUp.SimulatedVaf;
        }
        else
        {
            calcVaf = closestRatioDown.SimulatedVaf;
        }

        return calcVaf;
    }

    private class SimulatedVafCalcs
    {
        public final double DropoutRate;
        public final double SimulatedVaf;
        public final double ProbabilityTotalFrag1;
        public final double ProbabilityTotalFrag2Plus;

        public SimulatedVafCalcs(
                final double dropoutRate, final double simulatedVaf, final double probabilityTotalFrag1,
                final double probabilityTotalFrag2Plus)
        {
            DropoutRate = dropoutRate;
            SimulatedVaf = simulatedVaf;
            ProbabilityTotalFrag1 = probabilityTotalFrag1;
            ProbabilityTotalFrag2Plus = probabilityTotalFrag2Plus;
        }

        public double fragmentRatio() { return ProbabilityTotalFrag1 > 0 ? ProbabilityTotalFrag2Plus / ProbabilityTotalFrag1 : 0; }

        public String toString() { return format("simVaf(%.4f) doRate(%.2f) ratio(%.3f)", SimulatedVaf, DropoutRate, fragmentRatio()); }
    }

    public static BufferedWriter initialiseWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(DROPOUT_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(config.multipleSamples())
                writer.write("PatientId\t");

            writer.write("SampleId\tVariantCount\tSimulatedVaf\tDropoutRate");
            writer.write("\tProbTotalFrag1\tProbTotalFrag2Plus\tObsFrag1\tObsFrag2Plus");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeSimulatedDropoutData(
            final BufferedWriter writer, final PurityConfig config, final SampleData sample, final String sampleId,
            int varCount, final SimulatedVafCalcs simVafCalcs, int observedFrag1, int observedFrag2Plus)
    {
        if(writer == null)
            return;

        try
        {
            if(config.multipleSamples())
                writer.write(format("%s\t", sample.PatientId));

            writer.write(format("%s\t%d\t%.4f\t%.2f", sampleId, varCount, simVafCalcs.SimulatedVaf, simVafCalcs.DropoutRate));

            writer.write(format("\t%.6f\t%.6f\t%d\t%d",
                    simVafCalcs.ProbabilityTotalFrag1, simVafCalcs.ProbabilityTotalFrag2Plus, observedFrag1, observedFrag2Plus));
            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write dropout calc file: {}", e.toString());
        }
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
        /*

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
    public void calculate(final String sampleId)
    {
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

     */


}
