package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.DROPOUT_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityResult.INVALID_RESULT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.common.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class LowCountModel extends ClonalityModel
{
    public LowCountModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        super(config, resultsWriter, sample,  variants);
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
                else if(sampleFragData.AlleleCount >= 2)
                    ++observedFrag2Plus;
            }
        }

        if(observedFrag1 == 0)
            return INVALID_RESULT;

        List<SimulatedVafCalcs> simulatedVafCalcs = Lists.newArrayList();

        // now test each simluated dropout rate and VAF
        for(double dropoutRate = PurityConstants.DROPOUT_RATE_INCREMENT; dropoutRate < 0.95; dropoutRate += PurityConstants.DROPOUT_RATE_INCREMENT)
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
        double[] calcValues = findVafRatio(observedRatio, simulatedVafCalcs);
        double calcVaf = calcValues[0];
        double calcDropout = calcValues[1];

        double lowObservedRatio = max(observedFrag2Plus - 2, 0) / (double)observedFrag1;
        double lowCalcVaf = findVafRatio(lowObservedRatio, simulatedVafCalcs)[0];

        double highObservedRatio = (observedFrag2Plus + 2) / (double)observedFrag1;
        double highCalcVaf = findVafRatio(highObservedRatio, simulatedVafCalcs)[0];

        CT_LOGGER.debug(format("sample(%s) low-count model: obsRatio(%.2f 1=%d 2+=%d) vaf((%.6f low=%.6f high=%.6f)",
                sampleId, observedRatio, observedFrag1, observedFrag2Plus, calcVaf, lowCalcVaf, highCalcVaf));

        return new ClonalityResult(
                ClonalityMethod.LOW_COUNT, calcVaf, lowCalcVaf, highCalcVaf, observedFrag1 + observedFrag2Plus, calcDropout);
    }

    private static double[] findVafRatio(double observedRatio, final List<SimulatedVafCalcs> simulatedVafCalcs)
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
        double calcDropout = 0;

        if(closestRatioUp != null && closestRatioDown != null)
        {
            double upperFraction = (observedRatio - closestRatioDown.fragmentRatio())
                    / (closestRatioUp.fragmentRatio() - closestRatioDown.fragmentRatio());

            calcVaf = upperFraction * closestRatioUp.SimulatedVaf + (1 - upperFraction) * closestRatioDown.SimulatedVaf;
            calcDropout = upperFraction * closestRatioUp.DropoutRate + (1 - upperFraction) * closestRatioDown.DropoutRate;
        }
        else if(closestRatioUp != null)
        {
            calcVaf = closestRatioUp.SimulatedVaf;
            calcDropout = closestRatioUp.DropoutRate;
        }
        else
        {
            calcVaf = closestRatioDown.SimulatedVaf;
            calcDropout = closestRatioDown.DropoutRate;
        }

        return new double[] { calcVaf, calcDropout };
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
}
