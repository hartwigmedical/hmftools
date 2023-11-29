package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class CigarFrequency
{
    public final SAMRecord SampleRead;
    public int Frequency;

    public CigarFrequency(final SAMRecord sampleRead)
    {
        SampleRead = sampleRead;
        Frequency = 1;
    }

    public String toString()
    {
        return format("count(%d)", Frequency);
    }

    public static Map<String,CigarFrequency> buildFrequencies(final List<SAMRecord> reads)
    {
        return buildFrequencies(reads, false);
    }

    public static Map<String,CigarFrequency> buildMateFrequencies(final List<SAMRecord> reads)
    {
        return buildFrequencies(reads, true);
    }

    private static Map<String,CigarFrequency> buildFrequencies(final List<SAMRecord> reads, boolean useMateCigar)
    {
        Map<String,CigarFrequency> cigarFrequencies = Maps.newHashMap();

        for(SAMRecord read : reads)
        {
            final String cigarStr = useMateCigar ? read.getStringAttribute(MATE_CIGAR_ATTRIBUTE) : read.getCigarString();

            if(useMateCigar && (cigarStr == null || cigarStr.isEmpty()))
                continue;

            CigarFrequency frequency = cigarFrequencies.get(cigarStr);

            if(frequency == null)
                cigarFrequencies.put(cigarStr, new CigarFrequency(read));
            else
                ++frequency.Frequency;
        }

        return cigarFrequencies;
    }
}
