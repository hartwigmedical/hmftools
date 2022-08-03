package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.svprep.reads.ReadType.NO_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.rank;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Chromosome;
    public final int[] Positions;

    public String MateChromosome;
    public int MatePosStart;

    private final SAMRecord mRecord;
    private int mFragmentInsertSize;
    private final SupplementaryReadData mSupplementaryAlignment;

    private int mFilters;
    private ReadType mReadType;
    private boolean mWritten;

    public static ReadRecord from(final SAMRecord record) { return new ReadRecord(record); }

    public static final String UNMAPPED_CHR = "-1";

    public ReadRecord(final SAMRecord record)
    {
        mRecord = record;

        if(!record.getReadUnmappedFlag())
        {
            Chromosome = record.getReferenceName();
            Positions = new int[] { record.getStart(), record.getEnd() };
        }
        else
        {
            Chromosome = UNMAPPED_CHR;
            Positions = new int[] { 0, 0 };
        }

        if(!record.getMateUnmappedFlag())
        {
            MateChromosome = record.getMateReferenceName();
            MatePosStart = record.getMateAlignmentStart();
        }
        else
        {
            MateChromosome = UNMAPPED_CHR;
            MatePosStart = 0;
        }

        mFragmentInsertSize = abs(record.getInferredInsertSize());
        mSupplementaryAlignment = SupplementaryReadData.from(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
        mFilters = 0;
        mReadType = NO_SUPPORT;
        mWritten = false;
    }

    public String id() { return mRecord.getReadName(); }
    public final SAMRecord record() { return mRecord; }
    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }

    public byte orientation()
    {
        // first in pair has orientation of +1 if not reversed, and vice versa for the second in the pair
        if(isFirstOfPair())
            return !isReadReversed() ? POS_ORIENT : NEG_ORIENT;
        else
            return isReadReversed() ? NEG_ORIENT : POS_ORIENT;
    }

    public int flags() { return mRecord.getFlags(); }
    public Cigar cigar() { return mRecord.getCigar(); }
    public boolean isReadReversed() { return ( mRecord.getFlags() & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mRecord.getFlags() & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mRecord.getFlags() & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isUnmapped() { return (mRecord.getFlags() & SAMFlag.READ_UNMAPPED.intValue()) != 0; }
    public boolean isMateUnmapped() { return (mRecord.getFlags() & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }

    public boolean hasFlag(final SAMFlag flag) { return (mRecord.getFlags() & flag.intValue()) != 0; }

    public SupplementaryReadData supplementaryAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null && HumanChromosome.contains(mSupplementaryAlignment.Chromosome); }

    public String readBases() { return mRecord.getReadString(); }
    public byte[] baseQualities() { return mRecord.getBaseQualities(); }

    public void setFilters(int filters) { mFilters = filters; }
    public int filters() { return mFilters; }

    public void setReadType(ReadType type) { setReadType(type, false); }

    public void setReadType(ReadType type, boolean checkRank)
    {
        if(!checkRank || rank(type) > rank(mReadType)) // keep the highest
            mReadType = type;
    }

    public ReadType readType() { return mReadType; }

    public void setWritten() { mWritten = true; }
    public boolean written() { return mWritten; }

    public short mapQuality() { return (short)mRecord.getMappingQuality(); }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public String toString()
    {
        return format("coords(%s:%d-%d) cigar(%s) mate(%s:%d) id(%s) flags(first=%s supp=%s reversed=%s) hasSupp(%s) type(%s)",
                Chromosome, start(), end(), cigar().toString(), MateChromosome, MatePosStart, id(),
                isFirstOfPair(), isSupplementaryAlignment(), isReadReversed(), mSupplementaryAlignment != null, mReadType);
    }

    public static int maxIndelLength(final Cigar cigar)
    {
        return cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.D || x.getOperator() == CigarOperator.I)
                .mapToInt(x -> x.getLength()).max().orElse(0);
    }

    public static String getSoftClippedBases(final SAMRecord record, final boolean isClippedLeft)
    {
        int scLength = isClippedLeft ? record.getCigar().getFirstCigarElement().getLength() : record.getCigar().getLastCigarElement().getLength();
        int readLength = record.getReadBases().length;
        int scStart = isClippedLeft ? 0 : readLength - scLength;
        int scEnd = isClippedLeft ? scLength : readLength;
        return record.getReadString().substring(scStart, scEnd);
    }

    public static boolean hasPolyATSoftClip(final ReadRecord read, final boolean isClippedLeft)
    {
        byte orientation = isClippedLeft ? NEG_ORIENT : POS_ORIENT;
        String scBases = getSoftClippedBases(read.record(), isClippedLeft);
        return isMobileLineElement(orientation, scBases);
    }

    public static int[] findIndelCoords(final ReadRecord read, int minIndelLength)
    {
        int maxIndelLength = maxIndelLength(read.cigar());

        if(maxIndelLength < minIndelLength)
            return null;

        // find the location of the internal delete or insert
        int indelStartPos = read.start() - 1;
        int indelEndPos = 0;
        for(CigarElement element : read.cigar())
        {
            if(element.getOperator() == M)
            {
                indelStartPos += element.getLength();
            }
            else if(element.getOperator() == D)
            {
                if(element.getLength() == maxIndelLength)
                {
                    indelEndPos = indelStartPos + element.getLength() + 1;
                    break;
                }

                indelStartPos += element.getLength();
            }
            else if(element.getOperator() == I)
            {
                if(element.getLength() == maxIndelLength)
                {
                    indelEndPos = indelStartPos + 1;
                    break;
                }
            }
            else
            {
                continue;
            }
        }

        if(indelEndPos <= indelStartPos)
            return null;

        return new int[] { indelStartPos, indelEndPos };
    }
}
