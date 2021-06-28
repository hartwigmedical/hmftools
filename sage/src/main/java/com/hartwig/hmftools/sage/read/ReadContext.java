package com.hartwig.hmftools.sage.read;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext
{
    public final int Position;
    public final String Repeat;
    public final int RepeatCount;
    public final String Microhomology;
    public final IndexedBases ReadBases;
    public final IndexedBases RefBases;

    private final boolean mIncompleteCore;

    private static final int BONUS_FLANK = 50;

    @VisibleForTesting
    ReadContext(final String repeat, final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, final byte[] readBases, final String microhomology)
    {

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length - 1);
        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        Position = refPosition;
        Repeat = repeat;
        Microhomology = microhomology;
        RepeatCount = 0;
        ReadBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases);
        RefBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases);
    }

    public ReadContext(final IndexedBases refBases, final IndexedBases readBases, final int repeatCount, final String repeat,
            final String microhomology)
    {
        if(refBases.Bases.length < readBases.Bases.length)
        {
            throw new IllegalArgumentException();
        }

        ReadBases = readBases;
        RefBases = refBases;
        Position = refBases.Position;
        mIncompleteCore = false;
        RepeatCount = repeatCount;
        Repeat = repeat;
        Microhomology = microhomology;
    }

    ReadContext(final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, @NotNull final IndexedBases refSequence,
            @NotNull final SAMRecord record)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, record.getReadBases().length - 1);

        Position = refPosition;
        Repeat = repeat;
        RepeatCount = repeatCount;
        Microhomology = microhomology;
        ReadBases = new IndexedBases(Position, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, record.getReadBases());

        int refIndex = refSequence.index(Position);
        RefBases = new IndexedBases(Position,
                refIndex,
                refIndex + adjLeftCentreIndex - readIndex,
                refIndex + adjRightCentreIndex - readIndex,
                0,
                refSequence.Bases);

        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(int leftCentreIndex, int rightCentreIndex, @NotNull final ReadContext readContext)
    {
        Position = readContext.Position;
        Repeat = readContext.Repeat;
        RepeatCount = readContext.RepeatCount;
        Microhomology = readContext.Microhomology;

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readContext.ReadBases.Bases.length - 1);
        int readIndex = readContext.ReadBases.Index;
        ReadBases = new IndexedBases(Position,
                readIndex,
                adjLeftCentreIndex,
                adjRightCentreIndex,
                readContext.ReadBases.FlankSize,
                readContext.readBases());

        int refIndex = readContext.RefBases.index(Position);
        RefBases = new IndexedBases(Position,
                refIndex,
                refIndex + adjLeftCentreIndex - readIndex,
                refIndex + adjRightCentreIndex - readIndex,
                0,
                readContext.RefBases.Bases);

        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(@NotNull final ReadContext clone)
    {
        Position = clone.Position;
        Repeat = clone.Repeat;
        RepeatCount = clone.RepeatCount;
        Microhomology = clone.Microhomology;
        mIncompleteCore = clone.mIncompleteCore;

        RefBases = IndexedBases.resize(Position,
                clone.RefBases.Index,
                clone.RefBases.LeftCoreIndex,
                clone.RefBases.RightCoreIndex,
                clone.RefBases.FlankSize,
                0,
                clone.RefBases.Bases);

        ReadBases = IndexedBases.resize(Position,
                clone.ReadBases.Index,
                clone.ReadBases.LeftCoreIndex,
                clone.ReadBases.RightCoreIndex,
                clone.ReadBases.FlankSize,
                BONUS_FLANK,
                clone.ReadBases.Bases);
    }

    @NotNull
    public ReadContext extend(int leftCentreIndex, int rightCentreIndex)
    {
        return new ReadContext(leftCentreIndex, rightCentreIndex, this);
    }

    @NotNull
    public ReadContext minimiseFootprint()
    {
        return new ReadContext(this);
    }

    public boolean incompleteFlanks()
    {
        return !ReadBases.flanksComplete();
    }

    public boolean incompleteCore()
    {
        return mIncompleteCore;
    }

    int avgCentreQuality(int readIndex, @NotNull final SAMRecord record)
    {
        int leftOffset = this.readIndex() - ReadBases.LeftCoreIndex;
        int rightOffset = ReadBases.RightCoreIndex - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        float quality = 0;
        for(int i = leftIndex; i <= rightIndex; i++)
        {
            quality += record.getBaseQualities()[i];
        }
        return Math.round(quality / (rightIndex - leftIndex + 1));
    }

    public boolean phased(int offset, @NotNull final ReadContext other)
    {
        return ReadBases.phased(offset, other.ReadBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        return ReadBases.isCentreCovered(otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(boolean wildcardAllowedInCoreMatch, int otherReadIndex, byte[] otherBases)
    {
        return ReadBases.matchAtPosition(wildcardAllowedInCoreMatch, otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(@NotNull final ReadContext other)
    {
        return ReadBases.matchAtPosition(false, other.ReadBases);
    }

    public int readBasesPositionIndex()
    {
        return ReadBases.Index;
    }

    public int readBasesLeftFlankIndex()
    {
        return ReadBases.LeftFlankIndex;
    }

    public int readBasesRightFlankIndex()
    {
        return ReadBases.RightFlankIndex;
    }

    public int readBasesLeftCentreIndex()
    {
        return ReadBases.LeftCoreIndex;
    }

    public int readBasesRightCentreIndex()
    {
        return ReadBases.RightCoreIndex;
    }

    @Override
    public String toString()
    {
        return ReadBases.centerString();
    }

    @VisibleForTesting
    @NotNull
    public String centerBases()
    {
        return ReadBases.centerString();
    }

    @NotNull
    public String microhomology()
    {
        return Microhomology;
    }

    public byte[] readBases()
    {
        return ReadBases.Bases;
    }

    private int readIndex()
    {
        return ReadBases.Index;
    }

    public int maxFlankLength()
    {
        return ReadBases.maxFlankLength();
    }

    public byte[] refTrinucleotideContext(int position)
    {
        return RefBases.trinucleotideContext(position);
    }

    public int length()
    {
        return readBasesRightFlankIndex() - readBasesLeftFlankIndex() + 1;
    }

    @NotNull
    public String leftFlankString()
    {
        return ReadBases.leftFlankString();
    }

    @NotNull
    public String rightFlankString()
    {
        return ReadBases.rightFlankString();
    }

}
