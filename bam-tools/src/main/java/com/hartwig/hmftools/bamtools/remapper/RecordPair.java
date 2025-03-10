package com.hartwig.hmftools.bamtools.remapper;

import java.util.Objects;

import com.hartwig.hmftools.common.codon.Nucleotides;

import htsjdk.samtools.SAMRecord;

public class RecordPair
{
    public final SAMRecord First;
    public final SAMRecord Second;

    public RecordPair(final SAMRecord aRecord, final SAMRecord anotherRecord)
    {
        if(aRecord.getFirstOfPairFlag())
        {
            First = aRecord;
            Second = anotherRecord;
        }
        else if(aRecord.getSecondOfPairFlag())
        {
            Second = aRecord;
            First = anotherRecord;
        }
        else
        {
            throw new IllegalArgumentException("Records do not form a proper pair as neither is first.");
        }
    }

    public BamReadData leftData()
    {
        return BamReadData.fromRecord(First);
    }

    public BamReadData rightData()
    {
        return BamReadData.fromRecord(Second);
    }

    public byte[] leftBasesForRealignment()
    {
        return bases(First);
    }

    public byte[] rightBasesForRealignment()
    {
        return bases(Second);
    }

    private static byte[] bases(SAMRecord record)
    {
        if(record.getReadNegativeStrandFlag())
        {
            return Nucleotides.reverseComplementBases(record.getReadBases());
        }
        return record.getReadBases();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final RecordPair that = (RecordPair) o;
        return Objects.equals(First, that.First) && Objects.equals(Second, that.Second);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(First, Second);
    }
}
