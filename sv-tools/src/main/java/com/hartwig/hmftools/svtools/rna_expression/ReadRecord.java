package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    private final SAMRecord mSamRecord;

    public final String Id;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;

    public final String ReadBases;
    public final int Length; // of bases
    public final Cigar Cigar;

    public ReadRecord(SAMRecord record)
    {
        mSamRecord = record;

        Id = record.getReadName();
        Chromosome = record.getReferenceName();
        PosStart = record.getStart();
        PosEnd = record.getEnd();
        ReadBases = record.getReadString();
        Length = ReadBases.length();
        Cigar = record.getCigar();
    }

    public ReadRecord(final String id, final String chromosome, long posStart, long posEnd, final String readBases, final Cigar cigar)
    {
        mSamRecord = null;
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Length = ReadBases.length();
        Cigar = cigar;
    }

    public final SAMRecord samRecord() { return mSamRecord; }

    public long range() { return PosEnd - PosStart; }

    public long upperReadStart()
    {
        // make use of Cigar info if present
        if(range() == Length)
            return PosStart;

        if(Cigar != null && Cigar.getLastCigarElement().getOperator() == CigarOperator.M)
            return PosEnd - Cigar.getLastCigarElement().getLength() + 1;

        return max(PosStart, PosEnd - Length + 1);
    }

    public long lowerReadEnd()
    {
        // make use of Cigar info if present
        if(range() == Length)
            return PosEnd;

        if(Cigar != null && Cigar.getFirstCigarElement().getOperator() == CigarOperator.M)
            return PosStart + Cigar.getFirstCigarElement().getLength() - 1;

        return min(PosEnd, PosStart + Length - 1);
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d, range=%d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, range(), Length, Cigar != null ? Cigar.toString() : "", Id);
    }

}
