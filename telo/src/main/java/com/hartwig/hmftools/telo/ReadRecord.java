package com.hartwig.hmftools.telo;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Id;
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;

    public final String ReadBases;
    public final String BaseQualityString;
    public final int Length; // of bases
    public final Cigar Cigar;

    private int mFlags;
    private String mMateChromosome;
    private int mMatePosStart;

    private int mFragmentInsertSize;
    private String mSupplementaryAlignment;
    private short mMapQuality;

    private boolean mHasTeloContent;
    private boolean mCompleteGroup;

    private static final String SUPPLEMENTARY_ATTRIBUTE = "SA";

    public static ReadRecord from(final SAMRecord record)
    {
        final String readId = record.isSecondaryAlignment() ? String.format("%s_%s",
                record.getReadName(), record.getAttribute("HI")) : record.getReadName();

        ReadRecord read = new ReadRecord(
                readId, record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getBaseQualityString(), record.getCigar(), record.getInferredInsertSize(), record.getFlags(),
                record.getMateReferenceName(), record.getMateAlignmentStart());

        read.setSuppAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
        read.setMapQuality((short)record.getMappingQuality());
        return read;
    }

    public ReadRecord(
            final String id, final String chromosome, int posStart, int posEnd, final String readBases, final String baseQualityString,
            final Cigar cigar, int insertSize, int flags, final String mateChromosome, int matePosStart)
    {
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Length = ReadBases.length();
        Cigar = cigar;
        BaseQualityString = baseQualityString;

        mFlags = flags;
        mMateChromosome = mateChromosome;
        mMatePosStart = matePosStart;
        mFragmentInsertSize = insertSize;
        mSupplementaryAlignment = null;
        mMapQuality = 0;
        mHasTeloContent = false;
        mCompleteGroup = false;
    }

    public boolean hasTeloContent() { return mHasTeloContent; }
    public void setTeloContent(boolean toggle) { mHasTeloContent = toggle; }

    public boolean completeGroup() { return mCompleteGroup; }
    public void markCompleteGroup() { mCompleteGroup = true; }

    public byte orientation()
    {
        // first in pair has orientation of +1 if not reversed, and vice versa for the second in the pair
        if(isFirstOfPair())
            return !isReadReversed() ? 1 : (byte)-1;
        else
            return isReadReversed() ? (byte)-1 : 1;
    }

    public int flags() { return mFlags; }
    public boolean isReadReversed() { return (mFlags & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mFlags & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isDuplicate() { return (mFlags & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean isUnmapped() { return (mFlags & SAMFlag.READ_UNMAPPED.intValue()) != 0; }
    public boolean isMateUnmapped() { return (mFlags & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }

    public void setSuppAlignment(final String suppAlign) { mSupplementaryAlignment = suppAlign; }
    public String getSuppAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null; }

    public void setMapQuality(short mapQuality) { mMapQuality = mapQuality; }
    public short mapQuality() { return mMapQuality; }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public String mateChromosome() { return mMateChromosome; }
    public int mateStartPosition() { return mMatePosStart; }

    public boolean matches(final ReadRecord other)
    {
        return Id.equals(other.Id) && Cigar.toString().equals(other.Cigar.toString()) && PosStart == other.PosStart && PosEnd == other.PosEnd;
    }

    public void setFlag(SAMFlag flag, boolean toggle)
    {
        if(toggle)
            mFlags |= flag.intValue();
        else
            mFlags &= ~flag.intValue();
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, Length, Cigar != null ? Cigar.toString() : "", Id);
    }

}
