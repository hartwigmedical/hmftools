package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.IndexedBases;

import htsjdk.samtools.SAMRecord;

/**
 * Represents a Read Context.
 * <p>
 * The core is incomplete if the centre indices are clipped to the read bases.
 */
public class ReadContext_
{
    public final int Position;
    public final String Repeat;
    public final int RepeatCount;
    public final String Microhomology;

    private IndexedBases mReadBases;
    private boolean mIncompleteCore;

    public ReadContext_(
            int position, final String repeat, int repeatCount, final String microhomology, final IndexedBases readBases,
            final boolean incompleteCore)
    {
        Position = position;
        RepeatCount = repeatCount;
        Repeat = repeat;
        Microhomology = microhomology;

        mReadBases = readBases;
        mIncompleteCore = incompleteCore;
    }

    /**
     * Clips centre indices based on readBases. The core is incomplete if the centre indices are clipped. Create IndexedBases based on these
     * clipped indices.
     */
    public static ReadContext_ fromReadRecord(
            final String microhomology, int repeatCount, final String repeat, int refPosition, int readIndex,
            int leftCentreIndex, int rightCentreIndex, int flankSize, final byte[] readBases)
    {
        int clippedLeftCentreIndex = max(0, leftCentreIndex);
        int clippedRightCentreIndex = min(rightCentreIndex, readBases.length - 1);
        boolean incompleteCore = leftCentreIndex != clippedLeftCentreIndex || rightCentreIndex != clippedRightCentreIndex;

        IndexedBases indexedBases = new IndexedBases(
                refPosition, readIndex, clippedLeftCentreIndex, clippedRightCentreIndex, flankSize, readBases);

        return new ReadContext_(refPosition, repeat, repeatCount, microhomology, indexedBases, incompleteCore);
    }

    /**
     * Changes the core. Clips the indices to the length of the bases, and updated mReadBases. The core is incomplete if the centre indices are
     * clipped.
     */
    public void extendCore(int leftCentreIndex, int rightCentreIndex)
    {
        int clippedLeftCentreIndex = max(0, leftCentreIndex);
        int clippedRightCentreIndex = min(rightCentreIndex, readBases().length - 1);

        mReadBases = new IndexedBases(mReadBases.Position, mReadBases.Index, clippedLeftCentreIndex, clippedRightCentreIndex, mReadBases.FlankSize, readBases());
        mIncompleteCore = leftCentreIndex != clippedLeftCentreIndex || rightCentreIndex != clippedRightCentreIndex;
    }

    /**
     * Map the core indices to around readIndex and compute the average base quality of the quality around this.
     * Returns this average base quality rounded to an int.
     */
    public int avgCentreQuality(int readIndex, final SAMRecord record)
    {
        int leftIndex = readBasesLeftCentreIndex() + readIndex - mReadBases.Index;
        int rightIndex = readBasesRightCentreIndex() + readIndex - mReadBases.Index;

        double quality = 0.0;
        for (int i = leftIndex; i <= rightIndex; ++i)
        {
            quality += record.getBaseQualities()[i];
        }

        return (int) Math.round(quality / (rightIndex - leftIndex + 1));
    }

    /**
     * Returns whether the length of the left and right flanks is not equal to the flank size.
     */
    public boolean hasIncompleteFlanks()
    {
        return !mReadBases.flanksComplete();
    }

    public boolean hasIncompleteCore() { return mIncompleteCore; }

    public int readBasesPositionIndex() { return mReadBases.Index; }
    public int readBasesLeftFlankIndex() { return mReadBases.LeftFlankIndex; }
    public int readBasesRightFlankIndex() { return mReadBases.RightFlankIndex; }
    public int readBasesLeftCentreIndex() { return mReadBases.LeftCoreIndex; }
    public int readBasesRightCentreIndex() { return mReadBases.RightCoreIndex; }

    public String microhomology() { return Microhomology; }

    public IndexedBases indexedBases() { return mReadBases; }
    public byte[] readBases() { return mReadBases.Bases; }

    public int maxFlankLength() { return mReadBases.maxFlankLength(); }

    @VisibleForTesting
    public String coreString() { return mReadBases.coreString(); }

    public String leftFlankString() { return mReadBases.leftFlankString(); }
    public String rightFlankString() { return mReadBases.rightFlankString(); }

    @Override
    public String toString() { return mReadBases.coreString(); }
}
