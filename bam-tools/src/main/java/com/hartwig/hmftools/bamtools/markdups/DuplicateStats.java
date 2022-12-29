package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class DuplicateStats
{
    public int ReadCount;
    public int InterPartitionUnclear;

    public Map<Integer,Integer> DuplicateFrequencies;

    public DuplicateStats()
    {
        ReadCount = 0;
        InterPartitionUnclear = 0;
        DuplicateFrequencies = Maps.newHashMap();
    }

    public void merge(final DuplicateStats other)
    {
        ReadCount += other.ReadCount;
        InterPartitionUnclear += other.InterPartitionUnclear;

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
        if(fragments.stream().noneMatch(x -> x.status() == FragmentStatus.DUPLICATE))
            return;

        for(int i = 0; i < fragments.size() - 1;)
        {
            Fragment first = fragments.get(i);

            if(!first.status().isDuplicate())
            {
                ++i;
                continue;
            }

            int duplicateCount = 1;

            int j = i + 1;
            for(; j < fragments.size(); ++j)
            {
                Fragment second = fragments.get(j);

                if(!second.status().isDuplicate())
                    break;

                if(first.coordinates()[SE_START] == second.coordinates()[SE_START]
                && first.coordinates()[SE_END] == second.coordinates()[SE_END])
                {
                    ++duplicateCount;
                }
                else
                {
                    break;
                }
            }

            if(duplicateCount > 1)
            {
                addFrequency(duplicateCount);
            }

            i = j;
        }
    }
}
