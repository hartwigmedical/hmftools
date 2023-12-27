package com.hartwig.hmftools.esvee.common;

import java.util.Set;

import com.hartwig.hmftools.esvee.models.Record;

public class SampleSupport
{
    private final String mSampleName;
    private final boolean mIsGermline;
    private final int mQuality;
    private final Set<Record> mSplitReads;
    private final Set<Record> mDiscordantReads;
    private final int mSplitReadFragmentCount;
    private final int mDiscordantPairFragmentCount;

    public SampleSupport(
            final String sampleName, final boolean isGermline,
            final int quality, final Set<Record> splitReads, final Set<Record> discordantReads)
    {
        mSampleName = sampleName;
        mIsGermline = isGermline;
        mQuality = quality;
        mSplitReads = splitReads;
        mDiscordantReads = discordantReads;

        mSplitReadFragmentCount = (int) splitReads.stream().map(Record::getName).distinct().count();
        mDiscordantPairFragmentCount = (int) discordantReads.stream().map(Record::getName).distinct().count();
    }

    public String sampleName()
    {
        return mSampleName;
    }

    public boolean isGermline() { return mIsGermline; }

    public int quality()
    {
        return mQuality;
    }

    public Set<Record> splitReads()
    {
        return mSplitReads;
    }
    public Set<Record> discordantReads()
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
