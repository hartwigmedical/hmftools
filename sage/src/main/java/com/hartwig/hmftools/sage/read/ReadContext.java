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
    private final boolean mIncompleteCore;

    private static final int BONUS_FLANK = 50;

    @VisibleForTesting
    public ReadContext(
            final String repeat, final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, final byte[] readBases, final String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length - 1);
        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        Position = refPosition;
        Repeat = repeat;
        Microhomology = microhomology;
        RepeatCount = 0;
        mReadBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases);
    }

    public ReadContext(
            final IndexedBases refBases, final IndexedBases readBases, final int repeatCount, final String repeat, final String microhomology)
    {
        if(refBases.Bases.length < readBases.Bases.length)
        {
            throw new IllegalArgumentException();
        }

        Position = refBases.Position;
        mIncompleteCore = false;
        RepeatCount = repeatCount;
        Repeat = repeat;
        Microhomology = microhomology;
        mReadBases = readBases;
    }

    ReadContext(
            final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, final IndexedBases refSequence,
            final SAMRecord record)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, record.getReadBases().length - 1);

        Position = refPosition;
        Repeat = repeat;
        RepeatCount = repeatCount;
        Microhomology = microhomology;

        mReadBases = new IndexedBases(Position, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, record.getReadBases());

        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(int leftCentreIndex, int rightCentreIndex, final ReadContext readContext)
    {
        Position = readContext.Position;
        Repeat = readContext.Repeat;
        RepeatCount = readContext.RepeatCount;
        Microhomology = readContext.Microhomology;

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readContext.mReadBases.Bases.length - 1);
        int readIndex = readContext.mReadBases.Index;
        mReadBases = new IndexedBases(Position,
                readIndex,
                adjLeftCentreIndex,
                adjRightCentreIndex,
                readContext.mReadBases.FlankSize,
                readContext.readBases());

        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(@NotNull final ReadContext clone)
    {
        Position = clone.Position;
        Repeat = clone.Repeat;
        RepeatCount = clone.RepeatCount;
        Microhomology = clone.Microhomology;
        mIncompleteCore = clone.mIncompleteCore;

        mReadBases = IndexedBases.resize(Position,
                clone.mReadBases.Index,
                clone.mReadBases.LeftCoreIndex,
                clone.mReadBases.RightCoreIndex,
                clone.mReadBases.FlankSize,
                BONUS_FLANK,
                clone.mReadBases.Bases);
    }

    public ReadContext extend(int leftCentreIndex, int rightCentreIndex)
    {
        return new ReadContext(leftCentreIndex, rightCentreIndex, this);
    }

    public ReadContext minimiseFootprint()
    {
        return new ReadContext(this);
    }

    public boolean incompleteFlanks()
    {
        return !mReadBases.flanksComplete();
    }

    public boolean incompleteCore()
    {
        return mIncompleteCore;
    }

    public int avgCentreQuality(int readIndex, @NotNull final SAMRecord record)
    {
        int leftOffset = this.readIndex() - mReadBases.LeftCoreIndex;
        int rightOffset = mReadBases.RightCoreIndex - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        float quality = 0;
        for(int i = leftIndex; i <= rightIndex; i++)
        {
            quality += record.getBaseQualities()[i];
        }
        return Math.round(quality / (rightIndex - leftIndex + 1));
    }

    public boolean phased(int offset, final ReadContext other)
    {
        return mReadBases.phased(offset, other.mReadBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        return mReadBases.isCentreCovered(otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(boolean wildcardAllowedInCoreMatch, int otherReadIndex, byte[] otherBases)
    {
        return mReadBases.matchAtPosition(wildcardAllowedInCoreMatch, otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(@NotNull final ReadContext other)
    {
        return mReadBases.matchAtPosition(false, other.mReadBases);
    }

    public int readBasesPositionIndex()
    {
        return mReadBases.Index;
    }

    public int readBasesLeftFlankIndex()
    {
        return mReadBases.LeftFlankIndex;
    }

    public int readBasesRightFlankIndex()
    {
        return mReadBases.RightFlankIndex;
    }

    public int readBasesLeftCentreIndex()
    {
        return mReadBases.LeftCoreIndex;
    }

    public int readBasesRightCentreIndex()
    {
        return mReadBases.RightCoreIndex;
    }

    @Override
    public String toString()
    {
        return mReadBases.centerString();
    }

    @VisibleForTesting
    @NotNull
    public String centerBases() { return mReadBases.centerString(); }

    @NotNull
    public String microhomology() { return Microhomology; }

    public byte[] readBases() { return mReadBases.Bases; }

    private int readIndex() { return mReadBases.Index; }

    public int maxFlankLength() { return mReadBases.maxFlankLength(); }

    public int length()
    {
        return readBasesRightFlankIndex() - readBasesLeftFlankIndex() + 1;
    }

    @NotNull
    public String leftFlankString() { return mReadBases.leftFlankString(); }

    @NotNull
    public String rightFlankString() { return mReadBases.rightFlankString(); }

}
