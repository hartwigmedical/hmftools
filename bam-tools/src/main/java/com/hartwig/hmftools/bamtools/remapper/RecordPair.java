package com.hartwig.hmftools.bamtools.remapper;

import java.util.Objects;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RecordPair
{
    public final @NotNull SAMRecord First;
    public final @NotNull SAMRecord Second;

    public RecordPair(@NotNull final SAMRecord aRecord, @NotNull final SAMRecord anotherRecord)
    {
        if(!aRecord.getReadName().equals(anotherRecord.getReadName()))
        {
            throw new IllegalArgumentException("Left read name does not match right read name.");
        }
        if(!SamRecordUtils.firstInPair(aRecord))
        {
            throw new IllegalArgumentException("Left read is not first in pair.");
        }
        if(!SamRecordUtils.secondInPair(anotherRecord))
        {
            throw new IllegalArgumentException("Right read is not second in pair.");
        }
        if(aRecord.getProperPairFlag() != anotherRecord.getProperPairFlag())
        {
            throw new IllegalArgumentException("Proper pair values not equal.");
        }
        if(aRecord.getFirstOfPairFlag())
        {
            if(!anotherRecord.getSecondOfPairFlag())
            {
                throw new IllegalArgumentException("Records do not form a proper pair.");
            }
            this.First = aRecord;
            this.Second = anotherRecord;
        }
        else if(aRecord.getSecondOfPairFlag())
        {
            if(!anotherRecord.getFirstOfPairFlag())
            {
                throw new IllegalArgumentException("Records do not form a proper pair.");
            }
            this.Second = aRecord;
            this.First = anotherRecord;
        }
        else
        {
            throw new IllegalArgumentException("Records do not form a proper pair as neither is first.");
        }
    }

    public RawFastaData leftData()
    {
        return RawFastaData.fromRecord(First);
    }

    public RawFastaData rightData()
    {
        return RawFastaData.fromRecord(Second);
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
