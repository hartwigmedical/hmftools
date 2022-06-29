package com.hartwig.hmftools.svtools.sv_prep;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.generateMappedCoords;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MULTI_MAP_QUALITY_THRESHOLD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Id;
    public final String Chromosome;
    public final int[] Positions;

    public final String ReadBases;
    public final htsjdk.samtools.Cigar Cigar;

    private int mFlags;
    public String MateChromosome;
    public int MatePosStart;

    private final List<int[]> mMappedCoords;
    private int mFragmentInsertSize;
    private final SupplementaryReadData mSupplementaryAlignment;
    private short mMapQuality;
    private byte[] mBaseQualities;

    private static final String SUPPLEMENTARY_ATTRIBUTE = "SA";
    private static final String SECONDARY_ATTRIBUTE = "HI";

    public static ReadRecord from(final SAMRecord record)
    {
        final String readId = record.isSecondaryAlignment() ? String.format("%s_%s",
                record.getReadName(), record.getAttribute(SECONDARY_ATTRIBUTE)) : record.getReadName();

        ReadRecord read = new ReadRecord(
                readId, record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getCigar(), record.getInferredInsertSize(), record.getFlags(),
                record.getMateReferenceName(), record.getMateAlignmentStart(), record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

        read.setMapQuality((short)record.getMappingQuality());
        return read;
    }

    public ReadRecord(
            final String id, final String chromosome, int posStart, int posEnd, final String readBases, final Cigar cigar,
            int insertSize, int flags, final String mateChromosome, int matePosStart, final String suppData)
    {
        Id = id;
        Chromosome = chromosome;
        Positions = new int[] { posStart, posEnd };
        ReadBases = readBases;
        Cigar = cigar;
        MateChromosome = mateChromosome;
        MatePosStart = matePosStart;

        mFlags = flags;

        List<int[]> mappedCoords = generateMappedCoords(Cigar, Positions[SE_START]);
        mMappedCoords = Lists.newArrayListWithCapacity(mappedCoords.size());
        mMappedCoords.addAll(mappedCoords);

        mFragmentInsertSize = insertSize;

        mSupplementaryAlignment = SupplementaryReadData.from(suppData);
        mMapQuality = 0;
        mBaseQualities = null;
    }

    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }
    public int range() { return Positions[SE_END] - Positions[SE_START]; }

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
    public boolean isMateNegStrand() { return (mFlags & SAMFlag.MATE_REVERSE_STRAND.intValue()) != 0; }
    public boolean isMateUnmapped() { return (mFlags & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }
    public boolean isInversion() { return isReadReversed() == isMateNegStrand(); }
    public boolean isProperPair() { return (mFlags & SAMFlag.PROPER_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mFlags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isSecondaryAlignment() { return (mFlags & SAMFlag.SECONDARY_ALIGNMENT.intValue()) != 0; }

    public void setFragmentInsertSize(int size) { mFragmentInsertSize = size; }

    public SupplementaryReadData supplementaryAlignment() { return mSupplementaryAlignment; }
    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null; }

    public void setMapQuality(short mapQuality) { mMapQuality = mapQuality; }
    public short mapQuality() { return mMapQuality; }

    public void setBaseQualities(final byte[] qualities) { mBaseQualities = qualities; }
    public byte[] baseQualities() { return mBaseQualities; }

    public boolean isMultiMapped() { return mMapQuality <= MULTI_MAP_QUALITY_THRESHOLD; }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

}
