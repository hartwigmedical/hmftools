package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_DEFAULT_ERROR_RATE;
import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_NOISE_RATE;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.SHORTENED;
import static com.hartwig.hmftools.sage.quality.MsiJitterCalcs.findApplicableParams;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.quality.MsiModelParams;

import org.apache.commons.math3.distribution.BinomialDistribution;

import htsjdk.samtools.SAMRecord;

public class JitterData
{
    private int mLengthened;
    private int mShortened;
    private double mQualBoost;
    private boolean mWithinNoise;

    public JitterData()
    {
        mLengthened = 0;
        mShortened = 0;
        mQualBoost = 0;
        mWithinNoise = false;
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
    public boolean isWithinNoise() { return mWithinNoise; }

    public String toString() { return format("short(%d) long(%d)", mShortened, mLengthened); }

    public static JitterMatch checkJitter(final VariantReadContext readContext, final SAMRecord record, int readVarIndex)
    {
        if(readContext.AllRepeats.isEmpty())
            return JitterMatch.NONE;

        final byte[] readBases = record.getReadBases();
        final byte[] readQuals = record.getBaseQualities();

        // try each repeat in turn, lengthening and shortening it
        for(RepeatInfo repeat : readContext.AllRepeats)
        {
            JitterMatch jitterMatch = checkRepeatJitterMatch(repeat, readContext, readVarIndex, readBases, readQuals);

            if(jitterMatch != NONE)
                return jitterMatch;
        }

        return JitterMatch.NONE;
    }

    private static JitterMatch checkRepeatJitterMatch(
            final RepeatInfo repeat, final VariantReadContext readContext, int readVarIndex, final byte[] readBases, final byte[] readQuals)
    {
        int repeatLength = repeat.repeatLength();

        for(int i = 0; i <= 1; ++i)
        {
            JitterMatch jitterType = (i == 0) ? SHORTENED : LENGTHENED;

            int flankReadIndexStart = readVarIndex - readContext.leftLength();
            int flankReadIndexEnd = readVarIndex + readContext.rightLength() - 1;

            if(jitterType == SHORTENED)
                flankReadIndexStart += repeatLength;
            else
                flankReadIndexStart -= repeatLength;

            if(flankReadIndexStart < 0 || flankReadIndexEnd >= readBases.length)
                return JitterMatch.NONE;

            int readContextIndex = 0;
            int readIndex = flankReadIndexStart;
            boolean allMatched = true;
            boolean indexAdjusted = false;

            for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
            {
                if(!indexAdjusted)
                {
                    if(jitterType == SHORTENED && readContextIndex == repeat.endIndex() - repeatLength)
                    {
                        indexAdjusted = true;
                        readIndex -= repeatLength;
                    }
                    else if(jitterType == LENGTHENED && readContextIndex == repeat.endIndex() + 1)
                    {
                        indexAdjusted = true;
                        readContextIndex -= repeatLength;
                    }
                }

                if(readIndex >= readBases.length || readContextIndex >= readContext.ReadBases.length)
                    break;

                if(readBases[readIndex] != readContext.ReadBases[readContextIndex])
                {
                    // mismatch cannot be in the core
                    if(readContextIndex >= readContext.CoreIndexStart && readContextIndex <= readContext.CoreIndexEnd)
                    {
                        allMatched = false;
                        break;
                    }

                    // and must be low-qual
                    if(readQuals[readIndex] >= MATCHING_BASE_QUALITY)
                    {
                        allMatched = false;
                        break;
                    }
                }
            }

            if(allMatched)
                return jitterType;
        }

        return JitterMatch.NONE;
    }

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

        if(noiseOutcome == JitterNoiseOutcome.FILTER_VARIANT)
        {
            mWithinNoise = true;
            return;
        }

        mWithinNoise = false;
        mQualBoost = 1;

        if(noiseOutcome == JitterNoiseOutcome.SHORTENED_NOISE || noiseOutcome == JitterNoiseOutcome.BOTH_NOISE)
            mQualBoost += mShortened / (double)fullSupport;

        if(noiseOutcome == JitterNoiseOutcome.LENGTHENED_NOISE || noiseOutcome == JitterNoiseOutcome.BOTH_NOISE)
            mQualBoost += mLengthened / (double)fullSupport;
    }

    private enum JitterNoiseOutcome
    {
        SHORTENED_NOISE,
        LENGTHENED_NOISE,
        BOTH_NOISE,
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
            if(isWithinNoise(fullSupport, shortened, shortenedErrorRate))
                return JitterNoiseOutcome.FILTER_VARIANT;
        }
        else if(fullSupport > shortened)
        {
            double errorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount, -1);

            if(isWithinNoise(shortened, fullSupport, errorRate))
                outcome = JitterNoiseOutcome.SHORTENED_NOISE;
        }

        if(lengthened > fullSupport)
        {
            if(isWithinNoise(fullSupport, lengthened, lengthenedErrorRate))
                return JitterNoiseOutcome.FILTER_VARIANT;
        }
        else if(fullSupport > lengthened)
        {
            double errorRate = msiJitterCalcs.getErrorRate(allParams, maxRepeat.Bases, repeatCount, 1);

            if(isWithinNoise(lengthened, fullSupport, errorRate))
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

            if(prob > MSI_JITTER_NOISE_RATE)
                return JitterNoiseOutcome.FILTER_VARIANT;
        }

        return outcome;
    }

    private boolean isWithinNoise(int fullSupport, int jitterCount, double errorRate)
    {
        // a low full count relative to the total will be classified as within noise
        double jitterRatio = fullSupport / (double)(fullSupport + jitterCount);
        if(jitterRatio < 2 * errorRate)
            return true;

        // also test a p-value of jitter vs the full support counts
        BinomialDistribution distribution = new BinomialDistribution(fullSupport + jitterCount, errorRate);

        double prob = 1 - distribution.cumulativeProbability(min(fullSupport, jitterCount) - 1);

        return prob > MSI_JITTER_NOISE_RATE;
    }


    @VisibleForTesting
    public void setValues(int shortened, int lengthened)
    {
        mShortened = shortened;
        mLengthened = lengthened;
        mQualBoost = 0;
        mWithinNoise = false;
    }
}
