package com.hartwig.hmftools.esvee.old;

import java.util.Set;

import com.hartwig.hmftools.esvee.read.Read;

public class SampleSupport
{
    private final String mSampleName;
    private final int mQuality;
    private final Set<Read> mSplitReads;
    private final Set<Read> mDiscordantReads;
    private final int mSplitReadFragmentCount;
    private final int mDiscordantPairFragmentCount;

    public SampleSupport(
            final String sampleName, final int quality, final Set<Read> splitReads, final Set<Read> discordantReads)
    {
        mSampleName = sampleName;
        mQuality = quality;
        mSplitReads = splitReads;
        mDiscordantReads = discordantReads;

        mSplitReadFragmentCount = (int) splitReads.stream().map(Read::getName).distinct().count();
        mDiscordantPairFragmentCount = (int) discordantReads.stream().map(Read::getName).distinct().count();
    }

    public String sampleName()
    {
        return mSampleName;
    }

    public int quality()
    {
        return mQuality;
    }

    public Set<Read> splitReads()
    {
        return mSplitReads;
    }
    public Set<Read> discordantReads()
    {
        return mDiscordantReads;
    }

    public int splitReadFragmentCount()
    {
        return mSplitReadFragmentCount;
    }
    public int discordantPairFragmentCount() { return mDiscordantPairFragmentCount; }

    public int totalSupportFragmentCount()
    {
        return mSplitReadFragmentCount + mDiscordantPairFragmentCount;
    }
}
