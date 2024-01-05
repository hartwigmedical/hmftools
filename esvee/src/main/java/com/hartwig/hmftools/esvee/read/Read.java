package com.hartwig.hmftools.esvee.read;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import java.util.List;

import com.hartwig.hmftools.esvee.old.Alignment;
import com.hartwig.hmftools.esvee.old.Sequence;
import com.hartwig.hmftools.esvee.old.SequenceDecomposer;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class Read implements Sequence
{
    private final SAMRecord mRecord;
    private Read mMateRead;

    private List<Alignment> mAlignment;

    private List<SequenceDecomposer.Node> mDecomposition;

    public Read(final SAMRecord record)
    {
        mRecord = record;
        mAlignment = null;
        mDecomposition = null;
        mMateRead = null;
    }

    public SAMRecord bamRecord() { return mRecord; }

    public void setMateRead(final Read mate) { mMateRead = mate; }
    public boolean hasMateSet() { return mMateRead != null; }
    public Read mateRead() { return mMateRead; }

    @Override
    public String getName() { return mRecord.getReadName(); }

    public String getBasesString() { return mRecord.getReadString(); }

    @Override
    public byte[] getBases() { return mRecord.getReadBases(); }

    @Override
    public byte[] getBaseQuality() { return mRecord.getBaseQualities(); }

    public int getLength() { return mRecord.getReadLength(); }
    public int insertSize() { return mRecord.getInferredInsertSize(); }

    public String getChromosome() { return mRecord.getReferenceName(); }

    public Cigar getCigar() { return mRecord.getCigar(); }

    public int getAlignmentStart() { return mRecord.getAlignmentStart(); }
    public int getAlignmentEnd() { return mRecord.getAlignmentEnd(); }

    public int getUnclippedStart()  { return mRecord.getUnclippedStart(); }
    public int getUnclippedEnd() { return mRecord.getUnclippedEnd(); }

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
        if(isMateUnmapped())
            return getAlignmentEnd();

        if(mMateRead != null)
            return mMateRead.getAlignmentEnd();

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
        if(refPosition < getAlignmentStart())
        {
            if(!allowExtrapolation)
                return INVALID_INDEX;

            int baseDiff = getAlignmentStart() - refPosition;
            int softClipBases = leftSoftClipLength(getCigar());
            return baseDiff <= softClipBases ? softClipBases - baseDiff : INVALID_INDEX;
        }
        else if(refPosition > getAlignmentEnd())
        {
            if(!allowExtrapolation)
                return INVALID_INDEX;

            int baseDiff = refPosition - getAlignmentEnd();
            int softClipBases = rightSoftClipLength(getCigar());
            return baseDiff <= softClipBases ? bamRecord().getReadLength() - (softClipBases - baseDiff) - 1 : INVALID_INDEX;
        }

        return mRecord.getReadPositionAtReferencePosition(refPosition) - 1;
    }

    public Object getAttribute(final String name) { return mRecord.getAttribute(name); }

    public List<SequenceDecomposer.Node> decompose()
    {
        return mDecomposition == null ? (mDecomposition = SequenceDecomposer.decompose(this)) : mDecomposition;
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

    /*

    // public boolean isMateOnTheLeft() { return negativeStrand(); }

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

    public int hashCode()
    {
        final int firstOfPair = !isSecondOfPair() ? 1 : 0;
        return mRecord.getReadName().hashCode() ^ firstOfPair;
    }

    private synchronized List<Alignment> buildAlignmentBlocks()
    {
        if(mAlignment != null)
            return mAlignment;

        final List<Alignment> alignment = new ArrayList<>();
        if(mRecord.getReadUnmappedFlag())
        {
            alignment.add(Alignment.unmapped(mRecord.getReadLength()));
            return mAlignment = alignment;
        }

        final int mapQ = mRecord.getMappingQuality();
        int referencePosition = mRecord.getAlignmentStart();
        int readPosition = 1;

        for(CigarElement element : mRecord.getCigar().getCigarElements())
        {
            switch(element.getOperator())
            {
                case M:
                    alignment.add(new Alignment(mRecord.getReferenceName(), referencePosition, readPosition, element.getLength(), false, mapQ));
                    break;
                case I:
                    alignment.add(new Alignment("*", 0, readPosition, element.getLength(), false, mapQ));
                    break;
                case S:
                    alignment.add(new Alignment("?", 0, readPosition, element.getLength(), false, mapQ));
                    break;
            }

            if(element.getOperator().consumesReadBases())
                readPosition += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                referencePosition += element.getLength();
        }

        return mAlignment = alignment;
    }

    public List<Alignment> getAlignmentBlocks()
    {
        if(mAlignment == null)
            return buildAlignmentBlocks();

        return mAlignment;
    }
    */
}
