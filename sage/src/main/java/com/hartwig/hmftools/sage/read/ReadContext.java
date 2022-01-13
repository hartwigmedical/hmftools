package com.hartwig.hmftools.sage.read;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext
{
    public final int Position;
    public final String Repeat;
    public final int RepeatCount;
    public final String Microhomology;

    private IndexedBases mReadBases;

    // base quals are recorded across the flanks and core
    private final int[] mBaseQualities;

    // distance in bases from read index to start of base-qualities, cached in case the flank index changes
    private final int mBaseQualIndexOffset;

    private boolean mIncompleteCore;

    private static final int BONUS_FLANK = 50;

    public ReadContext(
            int position, final String repeat, final int repeatCount, final String microhomology, final IndexedBases readBases,
            final int[] baseQualities, final boolean incompleteCore)
    {
        Position = position;
        RepeatCount = repeatCount;
        Repeat = repeat;
        Microhomology = microhomology;

        mReadBases = readBases;

        mBaseQualities = baseQualities;
        mBaseQualIndexOffset = mReadBases.Index - mReadBases.LeftFlankIndex;

        mIncompleteCore = incompleteCore;
    }

    public static ReadContext fromReadRecord(
            final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, final SAMRecord record)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, record.getReadBases().length - 1);

        IndexedBases readBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, record.getReadBases());

        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        // extract base qualities for the flanks and core
        int coreFlankLength = readBases.leftFlankString().length() + readBases.centerString().length() + readBases.rightFlankString().length();

        int[] baseQualities = new int[coreFlankLength];
        final byte[] baseQuals = record.getBaseQualities();
        int startIndex = readBases.LeftFlankIndex;

        for(int i = 0; i < baseQualities.length; ++i)
        {
            baseQualities[i] = baseQuals[startIndex + i];
        }

        return new ReadContext(refPosition, repeat, repeatCount, microhomology, readBases, baseQualities, incompleteCore);
    }

    public void extendCore(int leftCentreIndex, int rightCentreIndex)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, mReadBases.Bases.length - 1);

        IndexedBases newReadBases = new IndexedBases(Position,
                mReadBases.Index,
                adjLeftCentreIndex,
                adjRightCentreIndex,
                mReadBases.FlankSize,
                readBases());

        mReadBases = newReadBases;
        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    public ReadContext cloneAndExpand()
    {
        IndexedBases newReadBases = IndexedBases.resize(Position,
                mReadBases.Index,
                mReadBases.LeftCoreIndex,
                mReadBases.RightCoreIndex,
                mReadBases.FlankSize,
                BONUS_FLANK,
                mReadBases.Bases);

        return new ReadContext(Position, Repeat, RepeatCount, Microhomology, newReadBases, mBaseQualities, mIncompleteCore);
    }

    public int[] baseQualities() { return mBaseQualities; }

    public boolean hasIncompleteFlanks()
    {
        return !mReadBases.flanksComplete();
    }
    public boolean hasIncompleteCore() { return mIncompleteCore; }

    public int avgCentreQuality(int readIndex, final SAMRecord record)
    {
        int leftOffset = this.readIndex() - mReadBases.LeftCoreIndex;
        int rightOffset = mReadBases.RightCoreIndex - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        double quality = 0;
        int baseLength = rightIndex - leftIndex + 1;
        for(int i = leftIndex; i <= rightIndex; i++)
        {
            quality += record.getBaseQualities()[i];
        }

        return (int)Math.round(quality / baseLength);
    }

    public boolean phased(int offset, final ReadContext other)
    {
        return mReadBases.phased(offset, other.mReadBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        return mReadBases.isCentreCovered(otherReadIndex, otherBases);
    }

    public ReadContextMatch matchAtPosition(boolean wildcardAllowedInCoreMatch, int otherReadIndex, byte[] otherBases)
    {
        return mReadBases.matchAtPosition(wildcardAllowedInCoreMatch, otherReadIndex, otherBases);
    }

    public ReadContextMatch matchAtPosition(@NotNull final ReadContext other)
    {
        return mReadBases.matchAtPosition(false, other.mReadBases);
    }

    public int readBasesPositionIndex() { return mReadBases.Index; }
    public int readBasesLeftFlankIndex() { return mReadBases.LeftFlankIndex; }
    public int readBasesRightFlankIndex() { return mReadBases.RightFlankIndex; }
    public int readBasesLeftCentreIndex() { return mReadBases.LeftCoreIndex; }
    public int readBasesRightCentreIndex() { return mReadBases.RightCoreIndex; }

    public String microhomology() { return Microhomology; }

    public IndexedBases indexedBases() { return mReadBases; }
    public byte[] readBases() { return mReadBases.Bases; }
    private int readIndex() { return mReadBases.Index; }

    public int maxFlankLength() { return mReadBases.maxFlankLength(); }

    public int length()
    {
        return readBasesRightFlankIndex() - readBasesLeftFlankIndex() + 1;
    }

    @VisibleForTesting
    public String centerBases() { return mReadBases.centerString(); }

    public String leftFlankString() { return mReadBases.leftFlankString(); }
    public String rightFlankString() { return mReadBases.rightFlankString(); }

    public String toString() { return mReadBases.centerString(); }
}
