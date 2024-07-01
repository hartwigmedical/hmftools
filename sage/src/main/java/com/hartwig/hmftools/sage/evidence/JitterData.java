package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.JITTER_QUAL_BOOST_MAX_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_HARD_FILTER_NOISE_RATE;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_NOISE_RATE;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.quality.MsiModelParams;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class JitterData
{
    private int mLengthened;
    private int mShortened;
    private double mQualBoost;
    private boolean mFilterOnNoise;
    private boolean mHardFilterOnNoise;

    public JitterData()
    {
        mLengthened = 0;
        mShortened = 0;
        mQualBoost = 1;
        mFilterOnNoise = false;
        mHardFilterOnNoise = false;
    }

    public void update(final JitterMatch jitterMatch)
    {
        if(jitterMatch == JitterMatch.LENGTHENED)
            mLengthened++;
        else if(jitterMatch == JitterMatch.SHORTENED)
            mShortened++;
    }

    public int shortened() { return mShortened; }
    public int lengthened() { return mLengthened; }
    public int[] summary() { return new int[] { mShortened, mLengthened }; }

    public double qualBoost() { return mQualBoost; }
    public boolean hardFilterOnNoise() { return mHardFilterOnNoise; }
    public boolean filterOnNoise() { return mFilterOnNoise; }

    public String toString() { return format("short(%d) long(%d)", mShortened, mLengthened); }

    public void setJitterQualFilterState(final MsiJitterCalcs msiJitterCalcs, final ReadContextCounter readContextCounter)
    {
        if(readContextCounter.readContext().MaxRepeat == null)
            return;

        int fullSupport = readContextCounter.readSupportCounts().Full;

        JitterNoiseOutcome noiseOutcome = calcNoiseOutcome(
                msiJitterCalcs, readContextCounter.sampleId(), readContextCounter.readContext().MaxRepeat,
                fullSupport, mShortened, mLengthened);

        if(noiseOutcome == null)
            return;

        if(noiseOutcome == JitterNoiseOutcome.HARD_FILTER_VARIANT)
        {
            mFilterOnNoise = true;
            mHardFilterOnNoise = true;
            return;
        }
        else if(noiseOutcome == JitterNoiseOutcome.FILTER_VARIANT)
        {
            mFilterOnNoise = true;
            return;
        }

        mFilterOnNoise = false;
        mHardFilterOnNoise = false;
        mQualBoost = 1;

        if(noiseOutcome == JitterNoiseOutcome.SHORTENED_NOISE || noiseOutcome == JitterNoiseOutcome.BOTH_NOISE)
            mQualBoost += mShortened / (double)fullSupport;

        if(noiseOutcome == JitterNoiseOutcome.LENGTHENED_NOISE || noiseOutcome == JitterNoiseOutcome.BOTH_NOISE)
            mQualBoost += mLengthened / (double)fullSupport;

        mQualBoost = min(mQualBoost, JITTER_QUAL_BOOST_MAX_PERC);
    }

    private enum JitterNoiseOutcome
    {
        NOISE,
        SHORTENED_NOISE,
        LENGTHENED_NOISE,
        BOTH_NOISE,
        HARD_FILTER_VARIANT,
        FILTER_VARIANT;
    }

    private JitterNoiseOutcome calcNoiseOutcome(
            final MsiJitterCalcs msiJitterCalcs, final String sampleId, final RepeatInfo maxRepeat,
            int fullSupport, int shortened, int lengthened)
    {
        List<MsiModelParams> allParams = msiJitterCalcs.getSampleParams(sampleId);

        if(allParams == null)
            return null;

        int repeatCount = maxRepeat.Count;
        double shortenedErrorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount - 1, 1);
        double lengthenedErrorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount + 1, -1);

        JitterNoiseOutcome outcome = null;

        if(shortened > fullSupport)
        {
            JitterNoiseOutcome filterOutcome = calcNoiseOutcome(fullSupport, shortened, shortenedErrorRate);

            if(filterOutcome != JitterNoiseOutcome.NOISE)
                return filterOutcome;
        }
        else if(fullSupport > shortened)
        {
            double errorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount, -1);

            // note the values are swapped here - ie asking if full support looks like jitter vs the shortened jitter count
            if(calcNoiseOutcome(shortened, fullSupport, errorRate) != JitterNoiseOutcome.NOISE)
                outcome = JitterNoiseOutcome.SHORTENED_NOISE;
        }

        if(lengthened > fullSupport)
        {
            JitterNoiseOutcome filterOutcome = calcNoiseOutcome(fullSupport, lengthened, lengthenedErrorRate);

            if(filterOutcome != JitterNoiseOutcome.NOISE)
                return filterOutcome;
        }
        else if(fullSupport > lengthened)
        {
            double errorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount, 1);

            if(calcNoiseOutcome(lengthened, fullSupport, errorRate) != JitterNoiseOutcome.NOISE)
                outcome = outcome == JitterNoiseOutcome.SHORTENED_NOISE ? JitterNoiseOutcome.BOTH_NOISE : JitterNoiseOutcome.LENGTHENED_NOISE;
        }

        if(min(shortened, lengthened) >= fullSupport)
        {
            int total = fullSupport + shortened + lengthened;
            double jitterRatio = fullSupport / (double)total;
            double avgErrorRate = (shortenedErrorRate + lengthenedErrorRate) * 0.5;

            if(jitterRatio < 2 * avgErrorRate)
                return JitterNoiseOutcome.FILTER_VARIANT;

            BinomialDistribution distribution = new BinomialDistribution(total, avgErrorRate);

            double prob = 1 - distribution.cumulativeProbability(fullSupport - 1);

            if(prob > MSI_JITTER_HARD_FILTER_NOISE_RATE)
                return JitterNoiseOutcome.HARD_FILTER_VARIANT;
            else if(prob > MSI_JITTER_NOISE_RATE)
                return JitterNoiseOutcome.FILTER_VARIANT;
        }

        return outcome;
    }

    private JitterNoiseOutcome calcNoiseOutcome(int fullSupport, int jitterCount, double errorRate)
    {
        // checks whether the jitter count can be explained as noise vs the full count

        // test a p-value of jitter vs the full support counts
        BinomialDistribution distribution = new BinomialDistribution(fullSupport + jitterCount, errorRate);

        double prob = 1 - distribution.cumulativeProbability(min(fullSupport, jitterCount) - 1);

        if(prob > MSI_JITTER_HARD_FILTER_NOISE_RATE)
            return JitterNoiseOutcome.HARD_FILTER_VARIANT;
        else if(prob > MSI_JITTER_NOISE_RATE)
            return JitterNoiseOutcome.FILTER_VARIANT;

        // a low full count relative to the total will be classified as within noise
        double jitterRatio = fullSupport / (double)(fullSupport + jitterCount);
        if(jitterRatio < 2 * errorRate)
            return JitterNoiseOutcome.FILTER_VARIANT;

        return JitterNoiseOutcome.NOISE;
    }

    @VisibleForTesting
    public void setValues(int shortened, int lengthened)
    {
        mShortened = shortened;
        mLengthened = lengthened;
        mQualBoost = 1;
        mFilterOnNoise = false;
        mHardFilterOnNoise = false;
    }
}
