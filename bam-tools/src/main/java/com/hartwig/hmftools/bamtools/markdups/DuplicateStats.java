package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class DuplicateStats
{
    public int ReadCount;
    public int Duplicates;

    public Map<Integer,Integer> DuplicateFrequencies;

    public DuplicateStats()
    {
        ReadCount = 0;
        Duplicates = 0;
        DuplicateFrequencies = Maps.newHashMap();
    }

    public void merge(final DuplicateStats other)
    {
        ReadCount += other.ReadCount;
        Duplicates += other.Duplicates;

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

    public void addDuplicateFrequencies(final List<Fragment> fragments)
    {
        fragments.stream().filter(x -> x.status() == PRIMARY).forEach(x -> addFrequency(x.duplicateCount() + 1));
    }
}
