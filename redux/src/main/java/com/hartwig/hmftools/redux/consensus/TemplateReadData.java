package com.hartwig.hmftools.redux.consensus;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class TemplateReadData
{
    public final String ReadId;
    public final String Chromosome;
    public final int AlignmentStart;
    public final String Cigar;
    public final int Flags;

    public TemplateReadData(final String readId, final String chromosome, final int alignmentStart, final String cigar, final int flags)
    {
        ReadId = readId;
        Chromosome = chromosome;
        AlignmentStart = alignmentStart;
        Cigar = cigar;
        Flags = flags;
    }

    public boolean firstInPair() { return isFlagSet(SAMFlag.FIRST_OF_PAIR); }
    public boolean readNegativeStrandFlag() { return isFlagSet(SAMFlag.READ_REVERSE_STRAND); }
    public boolean mateNegativeStrandFlag() { return isFlagSet(SAMFlag.MATE_REVERSE_STRAND); }

    private boolean isFlagSet(final SAMFlag flag) { return (Flags & flag.intValue()) != 0; }

    public static TemplateReadData fromRead(final SAMRecord read)
    {
        return new TemplateReadData(
                read.getReadName(), read.getReferenceName(), read.getAlignmentStart(), read.getCigarString(), read.getFlags());
    }
}
