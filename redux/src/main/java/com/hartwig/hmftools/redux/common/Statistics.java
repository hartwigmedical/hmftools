package com.hartwig.hmftools.redux.common;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.redux.DuplicateFrequency;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.consensus.ConsensusStatistics;
import com.hartwig.hmftools.redux.duplicate.UmiStatistics;

public class Statistics
{
    // recorded by partition readers
    public long TotalReads;
    public long DuplicateReads;
    public long DuplicateGroups;

    public final Map<Integer, DuplicateFrequency> DuplicateFrequencies;

    public final UmiStatistics UmiStats;

    public final ConsensusStatistics ConsensusStats;

    public Statistics()
    {
        TotalReads = 0;
        DuplicateReads = 0;
        DuplicateGroups = 0;
        DuplicateFrequencies = Maps.newHashMap();
        UmiStats = new UmiStatistics();
        ConsensusStats = new ConsensusStatistics();
    }

    public void merge(final Statistics other)
    {
        TotalReads += other.TotalReads;
        DuplicateReads += other.DuplicateReads;
        DuplicateGroups += other.DuplicateGroups;

        for(DuplicateFrequency dupFreq : other.DuplicateFrequencies.values())
        {
            addFrequency(dupFreq.ReadCount, dupFreq.Count, dupFreq.DualStrandCount);
        }

        UmiStats.merge(other.UmiStats);
        ConsensusStats.merge(other.ConsensusStats);
    }

    public void addFrequency(int frequency)
    {
        addFrequency(frequency, 1, 0);
    }

    public void addNonDuplicateCounts(int fragmentCount)
    {
        if(fragmentCount > 0)
            addFrequency(1, fragmentCount, 0);
    }

    private void addFrequency(int duplicateCount, long count, long dualStrandCount)
    {
        int rounded = roundFrequency(duplicateCount);
        DuplicateFrequency dupFreq = DuplicateFrequencies.get(rounded);

        if(dupFreq == null)
        {
            dupFreq = new DuplicateFrequency(duplicateCount, 0, 0);
            DuplicateFrequencies.put(rounded, dupFreq);
        }

        dupFreq.Count += count;
        dupFreq.DualStrandCount += dualStrandCount;
    }

    public void addUmiGroup(final int fragmentCount, boolean hasDualStrand)
    {
        addFrequency(fragmentCount, 1, hasDualStrand ? 1 : 0);
        ++DuplicateGroups;
        DuplicateReads += fragmentCount - 1;
    }

    public void addDuplicateGroup(final int fragmentCount)
    {
        addFrequency(fragmentCount);
        ++DuplicateGroups;
        DuplicateReads += fragmentCount - 1; // excluding the primary fragment
    }

    public void logStats()
    {
        RD_LOGGER.info("stats: totalReads({}) duplicates({}) dupGroups({}) umiGroups({})",
                TotalReads, DuplicateReads, DuplicateGroups, UmiStats.UmiGroups);

        if(RD_LOGGER.isDebugEnabled())
        {
            RD_LOGGER.info("consensus stats: {}", ConsensusStats);

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            String dupFreqStr = frequencies.stream()
                    .map(x -> format("%d=%d", x, DuplicateFrequencies.get(x).Count))
                    .collect(Collectors.joining(", "));

            RD_LOGGER.debug("duplicate frequency: {}", dupFreqStr);
        }
    }

    public static int roundFrequency(int frequency)
    {
        if(frequency <= 10)
            return frequency;
        else if(frequency <= 100)
            return round(1.0f * frequency / 10) * 10;
        else if(frequency <= 1000)
            return round(1.0f * frequency / 100) * 100;
        else
            return round(1.0f * frequency / 1000) * 1000;
    }

    public void writeDuplicateStats(final ReduxConfig config)
    {
        try
        {
            String filename = DuplicateFrequency.generateFilename(config.OutputDir, config.SampleId);

            List<DuplicateFrequency> dupReadfrequencies = Lists.newArrayList();

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            for(Integer frequency : frequencies)
            {
                DuplicateFrequency dupFreq = DuplicateFrequencies.get(frequency);

                dupReadfrequencies.add(dupFreq);
            }

            DuplicateFrequency.write(filename, dupReadfrequencies, config.UMIs.Enabled);
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }
}
