package com.hartwig.hmftools.svtools.rna_expression;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    private final SAMRecord mSamRecord;

    public final String Id;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;
    public final int Length;

    public final String ReadBases;
    public final Cigar Cigar;

    public ReadRecord(SAMRecord record)
    {
        mSamRecord = record;

        Id = record.getReadName();
        Chromosome = record.getReferenceName();
        PosStart = record.getStart();
        PosEnd = record.getEnd();
        Length = record.getReadLength();
        ReadBases = record.getReadString();
        Cigar = record.getCigar();
    }

    public ReadRecord(final String id, final String chromosome, long posStart, long posEnd, final String readBases, final Cigar cigar)
    {
        mSamRecord = null;
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        Length = (int)(PosEnd - PosStart);
        ReadBases = readBases;
        Cigar = cigar;
    }

    public final SAMRecord samRecord() { return mSamRecord; }

    public long range() { return PosEnd - PosStart; }

    public String toString()
    {
        return String.format("%s: range(%s: %d -> %d, dist=%d) length(%d) cigar(%s)",
                Id, Chromosome, PosStart, PosEnd, PosEnd - PosStart, Length, Cigar != null ? Cigar.toString() : "");
    }

}
