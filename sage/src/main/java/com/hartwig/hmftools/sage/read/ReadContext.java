package com.hartwig.hmftools.sage.read;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext
{

    private static final int BONUS_FLANK = 50;

    private final int position;
    private final String repeat;
    private final int repeatCount;
    private final String microhomology;
    private final IndexedBases readBases;
    private final IndexedBases refBases;

    private final boolean incompleteCore;

    @VisibleForTesting
    ReadContext(final String repeat, final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, final byte[] readBases, final String microhomology)
    {

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length - 1);
        this.incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        this.position = refPosition;
        this.repeat = repeat;
        this.microhomology = microhomology;
        this.repeatCount = 0;
        this.readBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases);
        this.refBases = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases);
    }

    public ReadContext(final IndexedBases refBases, final IndexedBases readBases, final int repeatCount, final String repeat,
            final String microhomology)
    {
        if(refBases.bases().length < readBases.bases().length)
        {
            throw new IllegalArgumentException();
        }

        this.readBases = readBases;
        this.refBases = refBases;
        this.position = refBases.position();
        this.incompleteCore = false;
        this.repeatCount = repeatCount;
        this.repeat = repeat;
        this.microhomology = microhomology;
    }

    ReadContext(final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, @NotNull final IndexedBases refSequence,
            @NotNull final SAMRecord record)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, record.getReadBases().length - 1);

        this.position = refPosition;
        this.repeat = repeat;
        this.repeatCount = repeatCount;
        this.microhomology = microhomology;
        this.readBases = new IndexedBases(position, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, record.getReadBases());

        int refIndex = refSequence.index(position);
        this.refBases = new IndexedBases(position,
                refIndex,
                refIndex + adjLeftCentreIndex - readIndex,
                refIndex + adjRightCentreIndex - readIndex,
                0,
                refSequence.bases());

        this.incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(int leftCentreIndex, int rightCentreIndex, @NotNull final ReadContext readContext)
    {
        this.position = readContext.position;
        this.repeat = readContext.repeat;
        this.repeatCount = readContext.repeatCount;
        this.microhomology = readContext.microhomology;

        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readContext.readBases.bases().length - 1);
        int readIndex = readContext.readBases.index();
        this.readBases = new IndexedBases(position,
                readIndex,
                adjLeftCentreIndex,
                adjRightCentreIndex,
                readContext.readBases.flankSize(),
                readContext.readBases());

        int refIndex = readContext.refBases.index(position);
        this.refBases = new IndexedBases(position,
                refIndex,
                refIndex + adjLeftCentreIndex - readIndex,
                refIndex + adjRightCentreIndex - readIndex,
                0,
                readContext.refBases.bases());

        this.incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;
    }

    private ReadContext(@NotNull final ReadContext clone)
    {
        this.position = clone.position;
        this.repeat = clone.repeat;
        this.repeatCount = clone.repeatCount;
        this.microhomology = clone.microhomology;
        this.incompleteCore = clone.incompleteCore;

        this.refBases = IndexedBases.resize(position,
                clone.refBases.index(),
                clone.refBases.leftCentreIndex(),
                clone.refBases.rightCentreIndex(),
                clone.refBases.flankSize(),
                0,
                clone.refBases.bases());

        this.readBases = IndexedBases.resize(position,
                clone.readBases.index(),
                clone.readBases.leftCentreIndex(),
                clone.readBases.rightCentreIndex(),
                clone.readBases.flankSize(),
                BONUS_FLANK,
                clone.readBases.bases());
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
        return !readBases.flanksComplete();
    }

    public boolean incompleteCore()
    {
        return incompleteCore;
    }

    int minCentreQuality(int readIndex, SAMRecord record)
    {
        int leftOffset = this.readIndex() - readBases.leftCentreIndex();
        int rightOffset = readBases.rightCentreIndex() - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        int quality = Integer.MAX_VALUE;
        for(int i = leftIndex; i <= rightIndex; i++)
        {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    int avgCentreQuality(int readIndex, @NotNull final SAMRecord record)
    {
        int leftOffset = this.readIndex() - readBases.leftCentreIndex();
        int rightOffset = readBases.rightCentreIndex() - this.readIndex();

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
        return readBases.phased(offset, other.readBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        return readBases.isCentreCovered(otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(boolean wildcardAllowedInCoreMatch, int otherReadIndex, byte[] otherBases)
    {
        return readBases.matchAtPosition(wildcardAllowedInCoreMatch, otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(@NotNull final ReadContext other)
    {
        return readBases.matchAtPosition(false, other.readBases);
    }

    public int readBasesPositionIndex()
    {
        return readBases.index();
    }

    public int readBasesLeftFlankIndex()
    {
        return readBases.leftFlankIndex();
    }

    public int readBasesRightFlankIndex()
    {
        return readBases.rightFlankIndex();
    }

    public int readBasesLeftCentreIndex()
    {
        return readBases.leftCentreIndex();
    }

    public int readBasesRightCentreIndex()
    {
        return readBases.rightCentreIndex();
    }

    @Override
    public String toString()
    {
        return readBases.centerString();
    }

    @VisibleForTesting
    @NotNull
    public String centerBases()
    {
        return readBases.centerString();
    }

    @NotNull
    public String microhomology()
    {
        return microhomology;
    }

    @NotNull
    public String repeat()
    {
        return repeat;
    }

    public int repeatCount()
    {
        return repeatCount;
    }

    public byte[] readBases()
    {
        return readBases.bases();
    }

    private int readIndex()
    {
        return readBases.index();
    }

    private int flankSize()
    {
        return readBases.flankSize();
    }

    public int maxFlankLength()
    {
        return readBases.maxFlankLength();
    }

    public byte[] refTrinucleotideContext(int position)
    {
        return refBases.trinucleotideContext(position);
    }

    public int length()
    {
        return readBasesRightFlankIndex() - readBasesLeftFlankIndex() + 1;
    }

    public int coreLength()
    {
        return readBasesRightCentreIndex() - readBasesLeftCentreIndex() + 1;
    }

    @NotNull
    public String leftFlankString()
    {
        return readBases.leftFlankString();
    }

    @NotNull
    public String rightFlankString()
    {
        return readBases.rightFlankString();
    }

}
