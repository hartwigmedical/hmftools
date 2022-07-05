package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svprep.SvConstants.MULTI_MAP_QUALITY_THRESHOLD;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Chromosome;
    public final int[] Positions;

    public String MateChromosome;
    public int MatePosStart;
    public short MapQuality;

    private final SAMRecord mRecord;
    private int mFragmentInsertSize;
    private final SupplementaryReadData mSupplementaryAlignment;

    private int mFilters;

    private static final String SUPPLEMENTARY_ATTRIBUTE = "SA";
    private static final String SECONDARY_ATTRIBUTE = "HI";

    public static ReadRecord from(final SAMRecord record) { return new ReadRecord(record); }

    public ReadRecord(final SAMRecord record)
    {
        mRecord = record;
        Chromosome = record.getReferenceName();
        Positions = new int[] { record.getStart(), record.getEnd() };
        MateChromosome = record.getMateReferenceName();
        MatePosStart = record.getMateAlignmentStart();

        //final String readId = record.isSecondaryAlignment() ? format("%s_%s",
        //        record.getReadName(), record.getAttribute(SECONDARY_ATTRIBUTE)) : record.getReadName();

        mFragmentInsertSize = abs(record.getInferredInsertSize());
        mSupplementaryAlignment = SupplementaryReadData.from(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
        mFilters = 0;
    }

    public String id() { return mRecord.getReadName(); }
    public final SAMRecord record() { return mRecord; }
    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }

    public byte orientation()
    {
        // first in pair has orientation of +1 if not reversed, and vice versa for the second in the pair
        if(isFirstOfPair())
            return !isReadReversed() ? 1 : (byte)-1;
        else
            return isReadReversed() ? (byte)-1 : 1;
    }

    public int flags() { return mRecord.getFlags(); }
    public Cigar cigar() { return mRecord.getCigar(); }
    public boolean isReadReversed() { return ( mRecord.getFlags() & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mRecord.getFlags() & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mRecord.getFlags() & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }

    public boolean isDuplicate() { return (mRecord.getFlags() & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean hasFlag(final SAMFlag flag) { return (mRecord.getFlags() & flag.intValue()) != 0; }

    public SupplementaryReadData supplementaryAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null; }

    public String readBases() { return mRecord.getReadString(); }
    public byte[] baseQualities() { return mRecord.getBaseQualities(); }

    public void setFilters(int filters) { mFilters = filters; }
    public int filters() { return mFilters; }

    public short mapQuality() { return (short)mRecord.getMappingQuality(); }
    public boolean isMultiMapped() { return mapQuality() <= MULTI_MAP_QUALITY_THRESHOLD; }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public String toString()
    {
        return format("coords(%s:%d-%d) cigar(%s) mate(%s:%d) id(%s)",
                Chromosome, start(), end(), cigar().toString(), MateChromosome, MatePosStart, id());
    }

    public static int maxDeleteLength(final Cigar cigar)
    {
        return cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.D).mapToInt(x -> x.getLength()).max().orElse(0);
    }
}
