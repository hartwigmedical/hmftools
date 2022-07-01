package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svprep.SvConstants.MULTI_MAP_QUALITY_THRESHOLD;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Id;
    public final String Chromosome;
    public final int[] Positions;

    public String MateChromosome;
    public int MatePosStart;
    public short MapQuality;

    public final Cigar mCigar;
    private int mFlags;
    private int mFragmentInsertSize;
    private final SupplementaryReadData mSupplementaryAlignment;

    private String mReadBases;
    private byte[] mBaseQualities;

    private static final String SUPPLEMENTARY_ATTRIBUTE = "SA";
    private static final String SECONDARY_ATTRIBUTE = "HI";

    public static ReadRecord from(final SAMRecord record)
    {
        final String readId = record.isSecondaryAlignment() ? String.format("%s_%s",
                record.getReadName(), record.getAttribute(SECONDARY_ATTRIBUTE)) : record.getReadName();

        ReadRecord read = new ReadRecord(
                readId, record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getCigar(), record.getReadString(), record.getBaseQualities(), record.getInferredInsertSize(), record.getFlags(),
                record.getMateReferenceName(), record.getMateAlignmentStart(), record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE),
                (short)record.getMappingQuality());

        return read;
    }

    public ReadRecord(
            final String id, final String chromosome, int posStart, int posEnd, final Cigar cigar, final String readBases,
            final byte[] baseQualities,  int insertSize, int flags, final String mateChromosome, int matePosStart, final String suppData,
            short mapQual)
    {
        Id = id;
        Chromosome = chromosome;
        Positions = new int[] { posStart, posEnd };
        mReadBases = null;
        mCigar = cigar;
        MateChromosome = mateChromosome;
        MatePosStart = matePosStart;
        MapQuality = mapQual;

        mReadBases = readBases;
        mBaseQualities = baseQualities;
        mFlags = flags;
        mFragmentInsertSize = insertSize;
        mSupplementaryAlignment = SupplementaryReadData.from(suppData);
    }

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

    public int flags() { return mFlags; }
    public Cigar cigar() { return mCigar; }
    public boolean isReadReversed() { return (mFlags & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mFlags & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mFlags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }

    public boolean isDuplicate() { return (mFlags & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean hasFlag(final SAMFlag flag) { return (mFlags & flag.intValue()) != 0; }

    public SupplementaryReadData supplementaryAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null; }

    public String readBases() { return mReadBases; }
    public byte[] baseQualities() { return mBaseQualities; }

    public void clearBaseData()
    {
        mReadBases = "";
        mBaseQualities = null;
    }

    public boolean isMultiMapped() { return MapQuality <= MULTI_MAP_QUALITY_THRESHOLD; }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

}
