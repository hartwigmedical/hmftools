package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public class BlacklistRepeatMaskerRecord
{
    private final String mChromosome;
    private final int mPosStart;
    private final int mPosEnd;
    private final int mSampleCount;
    private final int mDepthMin;
    private final int mDepthMax;
    private final Optional<String> mRepeatType;
    private final Optional<String> mRepeatInfo;
    private final Optional<Integer> mRepeatPosStart;
    private final Optional<Integer> mRepeatPosEnd;
    private final Optional<Integer> mCount;
    private final Optional<String> mOtherInfo;

    public BlacklistRepeatMaskerRecord(@NotNull final String chromosome, final int posStart, final int posEnd, final int sampleCount,
            final int depthMin, final int depthMax, @NotNull final Optional<String> repeatType, @NotNull final Optional<String> repeatInfo,
            @NotNull final Optional<Integer> repeatPosStart, @NotNull final Optional<Integer> repeatPosEnd,
            @NotNull final Optional<Integer> count, @NotNull final Optional<String> otherInfo)
    {
        mChromosome = stripChrPrefix(chromosome);
        mPosStart = posStart;
        mPosEnd = posEnd;
        mSampleCount = sampleCount;
        mDepthMin = depthMin;
        mDepthMax = depthMax;
        mRepeatType = repeatType;
        mRepeatInfo = repeatInfo;
        mRepeatPosStart = repeatPosStart;
        mRepeatPosEnd = repeatPosEnd;
        mCount = count;
        mOtherInfo = otherInfo;

        validate();
    }

    private void validate()
    {
        int emptyCount = 0;

        if(mRepeatType.isEmpty())
        {
            emptyCount++;
        }
        if(mRepeatInfo.isEmpty())
        {
            emptyCount++;
        }
        if(mRepeatPosStart.isEmpty())
        {
            emptyCount++;
        }
        if(mRepeatPosEnd.isEmpty())
        {
            emptyCount++;
        }
        if(mCount.isEmpty())
        {
            emptyCount++;
        }
        if(mOtherInfo.isEmpty())
        {
            emptyCount++;
        }

        if(emptyCount > 0 && emptyCount < 6)
        {
            MD_LOGGER.error("RepeatType, RepeatInfo, RepeatPosStart, RepeatPosEnd, Count, and OtherInfo must be all empty or all non-empty: {}", toString());
            System.exit(1);
        }
    }

    public BlacklistRegion getBlacklistRegion()
    {
        return new BlacklistRegion(mChromosome, mPosStart, mPosEnd, mSampleCount, mDepthMin, mDepthMax);
    }

    public Optional<RepeatMasker> getRepeatMasker()
    {
        if(mRepeatType.isEmpty())
        {
            return Optional.empty();
        }

        return Optional.of(new RepeatMasker(mRepeatType.get(), mRepeatInfo.get(), mRepeatPosStart.get(), mRepeatPosEnd.get(), mCount.get(), mOtherInfo.get()));
    }

    public String getChromosome()
    {
        return mChromosome;
    }

    public int getPosStart()
    {
        return mPosStart;
    }

    public int getPosEnd()
    {
        return mPosEnd;
    }

    public int getSampleCount()
    {
        return mSampleCount;
    }

    public int getDepthMin()
    {
        return mDepthMin;
    }

    public int getDepthMax()
    {
        return mDepthMax;
    }

    public Optional<String> getRepeatType()
    {
        return mRepeatType;
    }

    public Optional<String> getRepeatInfo()
    {
        return mRepeatInfo;
    }

    public Optional<Integer> getRepeatPosStart()
    {
        return mRepeatPosStart;
    }

    public Optional<Integer> getRepeatPosEnd()
    {
        return mRepeatPosEnd;
    }

    public Optional<Integer> getCount()
    {
        return mCount;
    }

    public Optional<String> getOtherInfo()
    {
        return mOtherInfo;
    }

    @Override
    public String toString()
    {
        return "BlacklistRepeatMaskerRecord{" + "Chromosome='" + mChromosome + '\'' + ", PosStart=" + mPosStart + ", PosEnd=" + mPosEnd
                + ", SampleCount=" + mSampleCount + ", DepthMin=" + mDepthMin + ", DepthMax=" + mDepthMax + ", RepeatType=" + mRepeatType
                + ", RepeatInfo=" + mRepeatInfo + ", RepeatPosStart=" + mRepeatPosStart + ", RepeatPosEnd=" + mRepeatPosEnd + ", Count="
                + mCount + ", OtherInfo=" + mOtherInfo + '}';
    }
}
