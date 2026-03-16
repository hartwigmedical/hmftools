package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class PeakGnomadFrequenciesChecker
{
    private static final double GNOMAD_FREQUENCY_TOLERANCE = 0.15;
    private final CandidatePeak Peak;

    public PeakGnomadFrequenciesChecker(final CandidatePeak peak)
    {
        Peak = peak;
    }

    public boolean checkGnomadFrequencies(GnomadFrequencySupplier frequencySupplier, double expectedMean)
    {
        DescriptiveStatistics lowAStats = new DescriptiveStatistics();
        DescriptiveStatistics highAStats = new DescriptiveStatistics();
        for(PositionEvidence evidence : Peak.allCapturedPoints())
        {
            double gnomadFrequency = frequencySupplier.getFrequency(evidence.Chromosome, evidence.Position);
            if(evidence.vaf() < 0.5)
            {
                lowAStats.addValue(gnomadFrequency);
            }
            else
            {
                highAStats.addValue(gnomadFrequency);
            }
        }
        AMB_LOGGER.debug(format("Peak at: %.3f has low AF gnomad mean: %.3f,  high AF gnomad mean: %.3f", Peak.vaf(), lowAStats.getMean(), highAStats.getMean()));
        if(lowAStats.getN() == 0 || highAStats.getN() == 0)
        {
            return false;
        }
        if(Math.abs(lowAStats.getMean() - expectedMean) > GNOMAD_FREQUENCY_TOLERANCE)
        {
            return false;
        }
        if(Math.abs(highAStats.getMean() - expectedMean) > GNOMAD_FREQUENCY_TOLERANCE)
        {
            return false;
        }
        return true;
    }
}
