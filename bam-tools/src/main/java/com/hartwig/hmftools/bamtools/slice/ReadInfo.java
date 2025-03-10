package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadInfo
{
    public final boolean IsPrimary;
    public final boolean FirstInPair;
    public final String Contig;
    public final int AlignmentStart;

    public ReadInfo(final boolean isPrimary, final boolean firstInPair, final String contig, final int alignmentStart)
    {
        IsPrimary = isPrimary;
        FirstInPair = firstInPair;
        Contig = contig;
        AlignmentStart = alignmentStart;
    }

    public static ReadInfo fromRead(final SAMRecord read)
    {
        return new ReadInfo(
                !read.getSupplementaryAlignmentFlag(),
                read.getReadPairedFlag() ? read.getFirstOfPairFlag() : true,
                read.getReferenceName(), read.getAlignmentStart());
    }

    public static ReadInfo fromReadMate(final SAMRecord read)
    {
        int matePosition;
        String mateConfig;

        if(!read.getMateUnmappedFlag())
        {
            matePosition = read.getMateAlignmentStart();
            mateConfig = read.getMateReferenceName();
        }
        else if(read.getSupplementaryAlignmentFlag())
        {
            // the supplementary of a primary with an unmapped mate needs to take the mate position from its primary location
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            matePosition = suppData.Position;
            mateConfig = suppData.Chromosome;
        }
        else
        {
            matePosition = read.getAlignmentStart();
            mateConfig = read.getReferenceName();
        }

        return new ReadInfo(
                true,
                read.getReadPairedFlag() ? read.getSecondOfPairFlag() : true, // opposite of read, otherwise meaningless or unpaired
                mateConfig, matePosition);
    }

    public static ReadInfo fromSupplementaryData(final SupplementaryReadData suppData, boolean isPrimary, boolean firstInPair)
    {
        return new ReadInfo(isPrimary, firstInPair, suppData.Chromosome, suppData.Position);
    }

    public boolean matches(final ReadInfo other)
    {
        return IsPrimary == other.IsPrimary && FirstInPair == other.FirstInPair
                && AlignmentStart == other.AlignmentStart && Contig.equals(other.Contig);
    }

    public String toString()
    {
        return format("%s:%d %s %s", Contig, AlignmentStart, IsPrimary ? "primary" : "supp", FirstInPair ? "first" : "second");
    }
}
