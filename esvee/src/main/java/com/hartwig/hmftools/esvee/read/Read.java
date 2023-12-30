package com.hartwig.hmftools.esvee.read;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.esvee.sequence.Alignment;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.processor.SequenceDecomposer;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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

    /*
    @Override
    public List<Alignment> getAlignmentBlocks()
    {
        if(mAlignment == null)
            return buildAlignmentBlocks();

        return mAlignment;
    }
    */

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

    public boolean equals(final Object obj)
    {
        if(this == obj)
            return true;
        if(obj == null || obj.getClass() != Record.class)
            return false;

        final Record other = (Record) obj;
        if(!other.mRecord.getReadName().equals(mRecord.getReadName()))
            return false;

        return isFirstOfPair() == other.isFirstOfPair();
    }
    */

    /*
    public void setBases(final byte[] bases, final byte[] quals)
    {
        mRecord.setReadBases(bases);
        mRecord.setBaseQualities(quals);
        mDecomposition = null;
    }

    @Override
    public void setChromosome(final String chromosome)
    {
        mRecord.setReadUnmappedFlag(false);
        mRecord.setReferenceName(chromosome);
        mAlignment = null;
    }

    @Override
    public void setAlignmentStart(final int position)
    {
        mRecord.setReadUnmappedFlag(false);
        mRecord.setAlignmentStart(position);
        mAlignment = null;
    }

    @Override
    public void setCigar(final Cigar cigar)
    {
        mRecord.setCigar(cigar);
        mAlignment = null;
    }

    @Override
    public void setCigar(final String cigar)
    {
        mRecord.setCigarString(cigar);
        mAlignment = null;
    }

    @Override
    public void setPositiveStrand(final boolean isPositiveStrand)
    {
        mRecord.setReadNegativeStrandFlag(!isPositiveStrand);
    }
    */

    /*
    public boolean isGermline()
    {
        return "germline".equals(mRecord.getHeader().getAttribute("userTag"));
    }

    public String sampleName()
    {
        @Nullable
        final String readGroupName = mRecord.getStringAttribute("RG");
        if(readGroupName == null)
            return fallbackReadGroup("not finding RG tag");
        @Nullable
        final SAMReadGroupRecord readGroup = mRecord.getHeader().getReadGroup(readGroupName);
        if(readGroup == null)
            return fallbackReadGroup("not finding matching read-group for tag in file");
        return readGroup.getSample();
    }

    private String fallbackReadGroup(final String fallbackReason)
    {
        final List<SAMReadGroupRecord> readGroups = mRecord.getHeader().getReadGroups();
        if(readGroups.isEmpty())
            throw new IllegalStateException("Cannot determine Sample Name for " + getName() + ". No read groups in file");
        final String firstRGSample = readGroups.get(0).getSample();
        if(readGroups.stream().allMatch(rg -> rg.getSample().equals(firstRGSample)))
            return firstRGSample; // All read-groups in this file are from the same sample.

        throw new IllegalStateException("Cannot determine fallback Sample Name for file for " + getName() + " after " + fallbackReason);
    }
    */
}
