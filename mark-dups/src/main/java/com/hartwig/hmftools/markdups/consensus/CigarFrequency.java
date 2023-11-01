package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

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
        Map<String,CigarFrequency> cigarFrequencies = Maps.newHashMap();

        for(SAMRecord read : reads)
        {
            CigarFrequency frequency = cigarFrequencies.get(read.getCigarString());

            if(frequency == null)
                cigarFrequencies.put(read.getCigarString(), new CigarFrequency(read));
            else
                ++frequency.Frequency;
        }

        return cigarFrequencies;
    }
}
