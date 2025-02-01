package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class PrepRead
{
    public final String Chromosome;

    private int mAlignmentStart;
    private int mAlignmentEnd;
    private int mUnclippedStart;
    private int mUnclippedEnd;

    public String MateChromosome;
    public int MatePosStart;

    private final SAMRecord mRecord;
    private int mFragmentInsertSize;
    private final SupplementaryReadData mSupplementaryAlignment;

    // read filtering and evaluation state
    private final int mAlignedBaseLength;
    private final int mSoftClipLengthLeft;
    private final int mSoftClipLengthRight;
    private final int mMaxIndelLength;
    private boolean mCheckedIndelCoords;
    private IndelCoords mIndelCoords;
    private int mFilters;

    private ReadType mReadType; // junction classification

    private boolean mWritten; // a check to avoid a read being written again

    public static final String UNMAPPED_CHR = "-1";

    public PrepRead(final SAMRecord record)
    {
        mRecord = record;

        int alignedBaseLength = 0;
        int softClipLengthLeft = 0;
        int softClipLengthRight = 0;
        int maxIndelLength = 0;

        for(int i = 0; i < record.getCigar().getCigarElements().size(); ++i)
        {
            CigarElement element = record.getCigar().getCigarElements().get(i);

            switch(element.getOperator())
            {
                case M:
                    alignedBaseLength += element.getLength();
                    break;

                case S:
                    if(i == 0)
                        softClipLengthLeft = element.getLength();
                    else
                        softClipLengthRight = element.getLength();
                    break;

                case D:
                case I:
                    maxIndelLength = max(element.getLength(), maxIndelLength);
                    break;

                default:
                    break;
            }
        }

        mAlignedBaseLength = alignedBaseLength;
        mSoftClipLengthLeft = softClipLengthLeft;
        mSoftClipLengthRight = softClipLengthRight;
        mMaxIndelLength = maxIndelLength;

        if(!record.getReadUnmappedFlag())
        {
            Chromosome = record.getReferenceName();
            mAlignmentStart = record.getStart();
            mAlignmentEnd = record.getEnd();
            mUnclippedStart = mAlignmentStart - mSoftClipLengthLeft;
            mUnclippedEnd = mAlignmentEnd + mSoftClipLengthRight;
        }
        else
        {
            Chromosome = UNMAPPED_CHR;
            mAlignmentStart = 0;
            mAlignmentEnd = 0;
            mUnclippedStart = 0;
            mUnclippedEnd = 0;
        }

        if(!mateUnmapped(record) && record.getMateAlignmentStart() > 0)
        {
            MateChromosome = record.getMateReferenceName();
            MatePosStart = record.getMateAlignmentStart();
        }
        else
        {
            MateChromosome = UNMAPPED_CHR;
            MatePosStart = 0;
        }

        mFragmentInsertSize = inferredInsertSizeAbs(record);
        mSupplementaryAlignment = SupplementaryReadData.extractAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

        mCheckedIndelCoords = false;
        mIndelCoords = null;

        mFilters = 0;
        mReadType = ReadType.NO_SUPPORT;
        mWritten = false;
    }

    public static PrepRead from(final SAMRecord record) { return new PrepRead(record); }

    public String id() { return mRecord.getReadName(); }
    public final SAMRecord record() { return mRecord; }
    public int start() { return mAlignmentStart; }
    public int end() { return mAlignmentEnd; }

    public int unclippedStart()  { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }
    public boolean isLeftClipped() { return mSoftClipLengthLeft > 0; }
    public boolean isRightClipped() { return mSoftClipLengthRight > 0; }
    public int leftClipLength() { return mSoftClipLengthLeft; }
    public int rightClipLength() { return mSoftClipLengthRight; }
    public int maxIndelLength()  { return mMaxIndelLength; }
    public int alignedBaseLength()  { return mAlignedBaseLength; }

    public Orientation orientation() { return !isReadReversed() ? FORWARD : REVERSE; }
    public Orientation mateOrientation() { return !hasFlag(SAMFlag.MATE_REVERSE_STRAND) ? FORWARD : REVERSE; }

    public int flags() { return mRecord.getFlags(); }
    public Cigar cigar() { return mRecord.getCigar(); }
    public boolean isReadReversed() { return ( mRecord.getFlags() & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return firstInPair(mRecord); }
    public boolean isSupplementaryAlignment() { return (mRecord.getFlags() & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isUnmapped() { return (mRecord.getFlags() & SAMFlag.READ_UNMAPPED.intValue()) != 0; }

    public boolean hasMate() { return MatePosStart > 0; }
    public boolean isMateUnmapped() { return (mRecord.getFlags() & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }

    public boolean hasFlag(final SAMFlag flag) { return (mRecord.getFlags() & flag.intValue()) != 0; }

    public SupplementaryReadData supplementaryAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null && HumanChromosome.contains(mSupplementaryAlignment.Chromosome); }

    public String readBases() { return mRecord.getReadString(); }
    public byte[] baseQualities() { return mRecord.getBaseQualities(); }

    public void addFilter(final ReadFilterType filterType) { mFilters |= filterType.flag(); }
    public boolean hasFilter(final ReadFilterType filterType) { return (mFilters & filterType.flag()) != 0; }
    public void removefilter(final ReadFilterType filterType) { mFilters &= ~filterType.flag(); }
    public boolean unfiltered() { return mFilters == 0; }
    public int filters() { return mFilters; }

    public void setReadType(ReadType type) { setReadType(type, false); }

    public void setReadType(ReadType type, boolean checkRank)
    {
        if(!checkRank || ReadType.rank(type) > ReadType.rank(mReadType)) // keep the highest
            mReadType = type;
    }

    public ReadType readType() { return mReadType; }

    public void setWritten() { mWritten = true; }
    public boolean written() { return mWritten; }

    public short mapQuality() { return (short)mRecord.getMappingQuality(); }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public IndelCoords indelCoords()
    {
        if(!mCheckedIndelCoords)
        {
            mCheckedIndelCoords = true;

            if(mMaxIndelLength >= MIN_INDEL_SUPPORT_LENGTH)
                mIndelCoords = findIndelCoords(start(), cigar().getCigarElements(), MIN_INDEL_SUPPORT_LENGTH);
        }

        return mIndelCoords;
    }

    public String toString()
    {
        return format("coords(%s:%d-%d) cigar(%s) mate(%s:%d) id(%s) flags(first=%s supp=%s reversed=%s) hasSupp(%s) type(%s)",
                Chromosome, start(), end(), cigar().toString(), MateChromosome, MatePosStart, id(),
                isFirstOfPair(), isSupplementaryAlignment(), isReadReversed(), mSupplementaryAlignment != null, mReadType);
    }
}
