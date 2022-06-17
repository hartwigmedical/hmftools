package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMRecord;

public class ReadContext
{
    public final int Position;
    public final String Repeat;
    public final int RepeatCount;
    public final String Microhomology;

    private IndexedBases mReadBases;

    private boolean mIncompleteCore;

    public ReadContext(
            int position, final String repeat, final int repeatCount, final String microhomology, final IndexedBases readBases,
            final boolean incompleteCore)
    {
        Position = position;
        RepeatCount = repeatCount;
        Repeat = repeat;
        Microhomology = microhomology;

        mReadBases = readBases;
        mIncompleteCore = incompleteCore;
    }

    public static ReadContext fromReadRecord(
            final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, final SAMRecord record)
    {
        int adjLeftCentreIndex = max(leftCentreIndex, 0);
        int adjRightCentreIndex = min(rightCentreIndex, record.getReadBases().length - 1);

        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBases = new IndexedBases(
                refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, record.getReadBases());

        return new ReadContext(refPosition, repeat, repeatCount, microhomology, readBases, incompleteCore);
    }

    public void extendCore(int leftCentreIndex, int rightCentreIndex)
    {
        int adjLeftCentreIndex = max(leftCentreIndex, 0);
        int adjRightCentreIndex = min(rightCentreIndex, mReadBases.Bases.length - 1);

        IndexedBases newReadBases = new IndexedBases(Position,
                mReadBases.Index,
                adjLeftCentreIndex,
                adjRightCentreIndex,
                mReadBases.FlankSize,
                readBases());

        mReadBases = newReadBases;
        mIncompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

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

    @VisibleForTesting
    public String coreString() { return mReadBases.coreString(); }

    public String leftFlankString() { return mReadBases.leftFlankString(); }
    public String rightFlankString() { return mReadBases.rightFlankString(); }

    public String toString() { return mReadBases.coreString(); }
}
