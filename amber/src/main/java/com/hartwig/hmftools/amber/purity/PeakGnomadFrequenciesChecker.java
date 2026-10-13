package com.hartwig.hmftools.amber.purity;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.GNOMAD_FREQUENCY_TOLERANCE;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class PeakGnomadFrequenciesChecker
{
    private final CandidatePeak mPeak;

    public PeakGnomadFrequenciesChecker(final CandidatePeak peak)
    {
        mPeak = peak;
    }

    public boolean checkGnomadFrequencies(GnomadFrequencySupplier frequencySupplier, double expectedMean)
    {
        DescriptiveStatistics lowAStats = new DescriptiveStatistics();
        DescriptiveStatistics highAStats = new DescriptiveStatistics();

        for(PositionEvidence evidence : mPeak.allCapturedPoints())
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

        AMB_LOGGER.trace(format("peak(%.3f) Gnomad lowMean(%.3f) highMean(%.3f)", mPeak.vaf(), lowAStats.getMean(), highAStats.getMean()));

        if(lowAStats.getN() == 0 || highAStats.getN() == 0)
            return false;

        if(abs(lowAStats.getMean() - expectedMean) > GNOMAD_FREQUENCY_TOLERANCE)
            return false;

        if(abs(highAStats.getMean() - expectedMean) > GNOMAD_FREQUENCY_TOLERANCE)
            return false;

        return true;
    }
}
