package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.DuplicateFrequency.roundFrequency;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.MarkDupsConfig;
import com.hartwig.hmftools.markdups.umi.UmiStatistics;

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

    public final Map<Integer,DuplicateFrequency> DuplicateFrequencies;

    public final UmiStatistics UmiStats;

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
        UmiStats = new UmiStatistics();
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

        for(DuplicateFrequency dupFreq : other.DuplicateFrequencies.values())
        {
            addFrequency(dupFreq.DuplicateCount, dupFreq.Frequency, dupFreq.DualStrandFrequency);
        }

        UmiStats.merge(other.UmiStats);
    }

    public void addFrequency(int frequency)
    {
        addFrequency(frequency, 1, 0);
    }

    private void addFrequency(int duplicateCount, int count, int dualStrandCount)
    {
        int rounded = roundFrequency(duplicateCount);
        DuplicateFrequency dupFreq = DuplicateFrequencies.get(rounded);

        if(dupFreq == null)
        {
            dupFreq = new DuplicateFrequency(duplicateCount);
            DuplicateFrequencies.put(rounded, dupFreq);
        }

        dupFreq.Frequency += count;
        dupFreq.DualStrandFrequency += dualStrandCount;
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
                MD_LOGGER.debug("duplicate frequency({}={})", frequency, DuplicateFrequencies.get(frequency).Frequency);
            }
        }
    }

    public void writeDuplicateStats(final MarkDupsConfig config)
    {
        try
        {
            String filename = config.formFilename("duplicate_freq");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("DuplicateReadCount\tFrequency");

            if(config.UMIs.Enabled)
                writer.write("\tDualStrandFrequency");

            writer.newLine();

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            for(Integer frequency : frequencies)
            {
                DuplicateFrequency dupFreq = DuplicateFrequencies.get(frequency);

                writer.write(format("%d\t%d", frequency, dupFreq.Frequency));

                if(config.UMIs.Enabled)
                    writer.write(format("\t%d", dupFreq.DualStrandFrequency));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }
}
