package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.SAMRecord;

public class BamReadData
{
    public final String ReadName;
    public final byte[] Bases;
    public final byte[] Qualities;

    public static BamReadData fromRecord(SAMRecord read)
    {
        if(read.getReadNegativeStrandFlag())
        {
            return new BamReadData(read.getReadName(), Nucleotides.reverseComplementBases(read.getReadBases()), Arrays.reverseArray(read.getBaseQualities()));
        }
        else
        {
            return new BamReadData(read.getReadName(), read.getReadBases(), read.getBaseQualities());
        }
    }

    public BamReadData(final String readName, final byte[] bases, final byte[] qualities)
    {
        ReadName = readName;
        Bases = bases;
        Qualities = qualities;
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
