package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import com.hartwig.hmftools.sage.quality.QualityConfig;

public class JitterData
{
    private double mPenalty;
    private int mLengthened;
    private int mShortened;

    public JitterData()
    {
        mPenalty = 0;
        mLengthened = 0;
        mShortened = 0;
    }

    public void update(final RealignedContext jitterRealign, final QualityConfig config)
    {
        if(jitterRealign.Type == LENGTHENED || jitterRealign.Type == SHORTENED)
        {
            mPenalty += jitterPenalty(config, jitterRealign.RepeatCount);

            if(jitterRealign.Type == LENGTHENED)
                mLengthened++;
            else
                mShortened++;
        }
    }

    private static double jitterPenalty(final QualityConfig config, int repeatCount)
    {
        return config.JitterPenalty * Math.max(0, repeatCount - config.JitterMinRepeatCount);
    }

    public int penalty() { return (int)mPenalty; }

    public int[] summary() { return new int[] { mShortened, mLengthened, penalty() }; }

    public String toString() { return format("short(%d) long(%d) penalty(%.0f)", mShortened, mLengthened, mPenalty); }

}
