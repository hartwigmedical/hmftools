package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.round;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.umi.UmiGroup;

public class Statistics
{
    public long TotalReads;
    public long DuplicateReads;
    public long DuplicateGroups;
    public long UmiGroups;

    // technical metrics
    public long LocalComplete;
    public long Incomplete;
    public long InterPartition;
    public long MissingMateCigar;

    public Map<Integer,Integer> DuplicateFrequencies;

    public Statistics()
    {
        TotalReads = 0;
        DuplicateReads = 0;
        DuplicateGroups = 0;
        UmiGroups = 0;
        InterPartition = 0;
        LocalComplete = 0;
        Incomplete = 0;
        MissingMateCigar = 0;
        DuplicateFrequencies = Maps.newHashMap();
    }

    public void merge(final Statistics other)
    {
        TotalReads += other.TotalReads;
        DuplicateReads += other.DuplicateReads;
        DuplicateGroups += other.DuplicateGroups;
        UmiGroups += other.UmiGroups;
        LocalComplete += other.LocalComplete;
        Incomplete += other.Incomplete;
        InterPartition += other.InterPartition;
        MissingMateCigar += other.MissingMateCigar;

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

    public void addUmiGroups(final List<Fragment> fragments, final List<UmiGroup> umiGroups)
    {
    }

    public void addDuplicateGroup(final List<Fragment> fragments)
    {
        addFrequency(fragments.size());
        DuplicateReads += fragments.size() - 1; // excluding the primary fragment
    }

    public void logStats()
    {
        MD_LOGGER.info("stats: totalReads({}) duplicates({}) duplicationGroups({}) umiGroups({})",
                TotalReads, DuplicateReads, DuplicateGroups, UmiGroups);

        if(MD_LOGGER.isDebugEnabled())
        {
            MD_LOGGER.debug("stats: fragments(complete={} incomplete={} interPartition={}) missingMateCigar({})",
                    LocalComplete, Incomplete, InterPartition, MissingMateCigar);

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            for(Integer frequency : frequencies)
            {
                MD_LOGGER.debug("duplicate frequency({}={})", frequency, DuplicateFrequencies.get(frequency));
            }
        }

    }
}
