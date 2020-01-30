package com.hartwig.hmftools.svtools.rna_expression;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    private final SAMRecord mSamRecord;

    public final String Id;
    public final long PosStart;
    public final long PosEnd;

    public final String ReadBases;
    public final Cigar Cigar;

    public ReadRecord(SAMRecord record)
    {
        mSamRecord = record;

        Id = record.getReadName();
        PosStart = record.getStart();
        PosEnd = record.getEnd();
        ReadBases = record.getReadString();
        Cigar = record.getCigar();
    }

    public ReadRecord(final String id, long posStart, long posEnd, final String readBases, final Cigar cigar)
    {
        mSamRecord = null;
        Id = id;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Cigar = cigar;
    }

    public final SAMRecord samRecord() { return mSamRecord; }
}
