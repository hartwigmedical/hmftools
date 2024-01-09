package com.hartwig.hmftools.esvee.read;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.CigarUtils.cigarElementsFromStr;
import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class Read
{
    private final SAMRecord mRecord;
    private Read mMateRead;

    // cached state and adjusted properties of the read
    private String mCigarString;
    private List<CigarElement> mCigarElements;

    private int mAlignmentStart;
    private int mAlignmentEnd;
    private int mUnclippedStart;
    private int mUnclippedEnd;
    private Integer mNumberOfEvents;
    private byte[] mBases;
    private byte[] mBaseQuals;

    public Read(final SAMRecord record)
    {
        mRecord = record;
        mMateRead = null;

        mCigarString = record.getCigarString();
        mCigarElements = cigarElementsFromStr(mCigarString);

        setBoundaries();
        mNumberOfEvents = 0;
        mBases = null;
        mBaseQuals = null;
    }

    private void setBoundaries()
    {
        mAlignmentStart = mRecord.getAlignmentStart();
        mUnclippedStart = mAlignmentStart;
        int currentPosition = mAlignmentStart;

        for(int i = 0; i < mCigarElements.size(); ++i)
        {
            CigarElement element = mCigarElements.get(i);

            if(i == 0 && element.getOperator() == S)
                mUnclippedStart -= element.getLength();

            if(element.getOperator().consumesReferenceBases())
                currentPosition += element.getLength();

            if(i == mCigarElements.size() - 1)
            {
                mAlignmentEnd = currentPosition - 1;

                mUnclippedEnd = element.getOperator() == S ? mAlignmentEnd + element.getLength() : mAlignmentEnd;
            }
        }
    }

    public SAMRecord bamRecord() { return mRecord; }

    public void setMateRead(final Read mate) { mMateRead = mate; }
    public boolean hasMateSet() { return mMateRead != null; }
    public Read mateRead() { return mMateRead; }

    public String getName() { return mRecord.getReadName(); }

    public String chromosome() { return mRecord.getReferenceName(); }

    public List<CigarElement> cigarElements() { return mCigarElements; }
    public String cigarString() { return mCigarString; }

    @Deprecated
    public Cigar getCigar() { return mRecord.getCigar(); }

    public int alignmentStart() { return mAlignmentStart; }
    public int alignmentEnd() { return mAlignmentEnd; }

    public int unclippedStart()  { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }

    // convenience
    public boolean isLeftClipped() { return mUnclippedStart < mAlignmentStart; }
    public boolean isRightClipped() { return mUnclippedEnd > mAlignmentEnd; }
    public int leftClipLength() { return mAlignmentStart - mUnclippedStart; }
    public int rightClipLength() { return mUnclippedEnd - mAlignmentEnd; }

    public byte[] getBases() { return mBases != null ? mBases : mRecord.getReadBases(); }
    public byte[] getBaseQuality() { return mBaseQuals != null ? mBaseQuals : mRecord.getBaseQualities(); }
    // public String getBasesString() { return mRecord.getReadString(); }
    public int basesLength() { return mBases != null ? mBases.length : mRecord.getReadBases().length; }
    public int insertSize() { return mRecord.getInferredInsertSize(); }

    // flags
    public boolean isUnmapped() { return mRecord.getReadUnmappedFlag(); }
    public boolean isPairedRead() { return mRecord.getReadPairedFlag(); }
    public boolean isFirstOfPair() { return mRecord.getReadPairedFlag() && mRecord.getFirstOfPairFlag(); }

    public boolean positiveStrand() { return !mRecord.getReadNegativeStrandFlag(); }
    public boolean negativeStrand() { return mRecord.getReadNegativeStrandFlag(); }
    public byte orientation() { return mRecord.getReadNegativeStrandFlag() ? NEG_ORIENT : POS_ORIENT; }
    public boolean firstInPair() { return mRecord.getFirstOfPairFlag(); }
    public boolean secondInPair() { return mRecord.getReadPairedFlag() && mRecord.getSecondOfPairFlag(); }

    public int mappingQuality() { return mRecord.getMappingQuality(); }

    public String mateChromosome() { return isMateMapped() ? mRecord.getMateReferenceName() : null; }
    public int mateAlignmentStart() { return mRecord.getMateAlignmentStart(); }

    public int mateAlignmentEnd()
    {
        // if used then should be calculated and cached
        if(isMateUnmapped())
            return alignmentEnd();

        if(mMateRead != null)
            return mMateRead.alignmentEnd();

        return getMateAlignmentEnd(mRecord);
    }

    public boolean isMateMapped() { return mRecord.getReadPairedFlag() && !mRecord.getMateUnmappedFlag(); }
    public boolean isMateUnmapped() { return mRecord.getReadPairedFlag() && mRecord.getMateUnmappedFlag(); }

    public boolean matePositiveStrand() { return !mRecord.getMateNegativeStrandFlag(); }
    public boolean mateNegativeStrand() { return mRecord.getMateNegativeStrandFlag(); }

    public static final int INVALID_INDEX = -1;

    public int getReadIndexAtReferencePosition(final int refPosition)
    {
        return getReadIndexAtReferencePosition(refPosition, false);
    }

    public int getReadIndexAtReferencePosition(final int refPosition, boolean allowExtrapolation)
    {
        // finds the read index given a reference position, and extrapolates outwards from alignments as required
        if(refPosition < mAlignmentStart)
        {
            if(!allowExtrapolation)
                return INVALID_INDEX;

            int baseDiff = mAlignmentStart - refPosition;
            int softClipBases = mAlignmentStart - mUnclippedStart;
            return baseDiff <= softClipBases ? softClipBases - baseDiff : INVALID_INDEX;
        }
        else if(refPosition > mAlignmentEnd)
        {
            if(!allowExtrapolation)
                return INVALID_INDEX;

            int baseDiff = refPosition - mAlignmentEnd;
            int softClipBases = mUnclippedEnd - mAlignmentEnd;
            return baseDiff <= softClipBases ? basesLength() - (softClipBases - baseDiff) - 1 : INVALID_INDEX;
        }

        return mRecord.getReadPositionAtReferencePosition(refPosition) - 1;
    }

    public Object getAttribute(final String name) { return mRecord.getAttribute(name); }

    public int numberOfEvents()
    {
        if(mNumberOfEvents != null)
            return mNumberOfEvents;

        Object numOfEvents = mRecord.getAttribute(NUM_MUTATONS_ATTRIBUTE);

        mNumberOfEvents = numOfEvents != null ? (int)numOfEvents : 0;
        return mNumberOfEvents;
    }

    public String toString()
    {
        return readToString(mRecord);
    }

    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public String sampleName() { return mRecord.getHeader().getAttribute(BAM_HEADER_SAMPLE_ID_TAG); }

    // public boolean isMateOnTheLeft() { return negativeStrand(); }

    /*
    public int impliedFragmentLength()
    {
        if(isMateMapped())
        {
            if(isMateOnTheLeft())
            {
                return getUnclippedEnd() - mRecord.getMateAlignmentStart();
            }
            else
            {
                final int mateEnd = mRecord.getMateAlignmentStart() + getLength();
                return mateEnd - getUnclippedStart();
            }
        }
        else
        {
            return getUnclippedEnd() - getUnclippedStart();
        }
    }

    */
}
