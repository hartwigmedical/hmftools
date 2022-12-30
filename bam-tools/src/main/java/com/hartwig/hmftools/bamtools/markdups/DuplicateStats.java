package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class DuplicateStats
{
    public long ReadCount;
    public long Duplicates;
    public long NoMateCigar;

    public Map<Integer,Integer> DuplicateFrequencies;

    public DuplicateStats()
    {
        ReadCount = 0;
        Duplicates = 0;
        NoMateCigar = 0;
        DuplicateFrequencies = Maps.newHashMap();
    }

    public void merge(final DuplicateStats other)
    {
        ReadCount += other.ReadCount;
        Duplicates += other.Duplicates;
        NoMateCigar += other.NoMateCigar;

        for(Map.Entry<Integer,Integer> entry : other.DuplicateFrequencies.entrySet())
        {
            Integer count = DuplicateFrequencies.get(entry.getKey());
            DuplicateFrequencies.put(entry.getKey(), count == null ? entry.getValue() : count + entry.getValue());
        }
    }

    public void addFrequency(int frequency)
    {
        int rounded;

        if(frequency <= 10)
            rounded = frequency;
        else if(frequency <= 100)
            rounded = round(frequency/10) * 10;
        else if(frequency <= 1000)
            rounded = round(frequency/100) * 100;
        else
            rounded = round(frequency/1000) * 1000;

        Integer count = DuplicateFrequencies.get(rounded);
        DuplicateFrequencies.put(rounded, count == null ? 1 : count + 1);
    }

    public void addDuplicateInfo(final List<Fragment> fragments)
    {
        for(Fragment fragment : fragments)
        {
            if(fragment.status() == PRIMARY)
                addFrequency(fragment.duplicateCount() + 1);

            if(fragment.status().isDuplicate())
                ++Duplicates;
        }
    }
}
