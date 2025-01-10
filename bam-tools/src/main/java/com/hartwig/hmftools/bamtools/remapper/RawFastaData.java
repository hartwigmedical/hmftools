package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RawFastaData
{
    @NotNull
    public final String ReadName;
    @NotNull
    public final byte[] Bases;
    @NotNull
    public final byte[] Qualities;

    public static RawFastaData fromRecord(SAMRecord read)
    {
        if(read.getReadNegativeStrandFlag())
        {
            return new RawFastaData(read.getReadName(), Nucleotides.reverseComplementBases(read.getReadBases()), Arrays.reverseArray(read.getBaseQualities()));
        }
        else
        {
            return new RawFastaData(read.getReadName(), read.getReadBases(), read.getBaseQualities());
        }
    }

    public RawFastaData(@NotNull final String readName, @NotNull final byte[] bases, @NotNull final byte[] qualities)
    {
        this.ReadName = readName;
        this.Bases = bases;
        this.Qualities = qualities;
    }

    @Override
    public String toString()
    {
        return "RawFastaData{" +
                "readName='" + ReadName + '\'' +
                ", bases=" + java.util.Arrays.toString(Bases) +
                ", qualities=" + java.util.Arrays.toString(Qualities) +
                '}';
    }
}
