package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.DuplicateFrequency.roundFrequency;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.consensus.ConsensusStatistics;
import com.hartwig.hmftools.redux.umi.UmiStatistics;

public class Statistics
{
    public long TotalReads;
    public long DuplicateReads;
    public long DuplicateGroups;

    // technical metrics
    public long LocalComplete; // fragments where all reads are in the same partition
    public long Incomplete; // read in same partition as base partition but not resolved immediately (eg an earlier supplementary)
    public long InterPartition; // reads where base partition (lower of mate or supplementary's primary) isn't the current partition
    public long MissingMateCigar;
    public long Unmapped; // fully, ie primary and mate
    public long PairedAltChromosome; // paired with a non-human chromosome

    public final Map<Integer,DuplicateFrequency> DuplicateFrequencies;

    public final UmiStatistics UmiStats;

    public final ConsensusStatistics ConsensusStats;

    public Statistics()
    {
        TotalReads = 0;
        DuplicateReads = 0;
        DuplicateGroups = 0;
        InterPartition = 0;
        LocalComplete = 0;
        Incomplete = 0;
        MissingMateCigar = 0;
        Unmapped = 0;
        PairedAltChromosome = 0;
        DuplicateFrequencies = Maps.newHashMap();
        UmiStats = new UmiStatistics();
        ConsensusStats = new ConsensusStatistics();
    }

    public void merge(final Statistics other)
    {
        TotalReads += other.TotalReads;
        DuplicateReads += other.DuplicateReads;
        DuplicateGroups += other.DuplicateGroups;
        LocalComplete += other.LocalComplete;
        Incomplete += other.Incomplete;
        InterPartition += other.InterPartition;
        MissingMateCigar += other.MissingMateCigar;
        Unmapped += other.Unmapped;
        PairedAltChromosome += other.PairedAltChromosome;

        for(DuplicateFrequency dupFreq : other.DuplicateFrequencies.values())
        {
            addFrequency(dupFreq.DuplicateCount, dupFreq.Frequency, dupFreq.DualStrandFrequency);
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

    private void addFrequency(int duplicateCount, long count, int dualStrandCount)
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
        RD_LOGGER.info("stats: totalReads({}) duplicates({}) duplicationGroups({}) umiGroups({}) {}",
                TotalReads, DuplicateReads, DuplicateGroups, UmiStats.UmiGroups, ConsensusStats);

        if(RD_LOGGER.isDebugEnabled())
        {
            RD_LOGGER.debug("stats: fragments(complete={} incomplete={} interPartition={} unmapped={} pairedAltChr={}))",
                    LocalComplete, Incomplete, InterPartition, Unmapped, PairedAltChromosome);

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            String dupFreqStr = frequencies.stream()
                    .map(x -> format("%d=%d", x, DuplicateFrequencies.get(x).Frequency))
                    .collect(Collectors.joining(", "));

            RD_LOGGER.debug("duplicate frequency: {}", dupFreqStr);
        }

        if(MissingMateCigar > 0)
        {
            RD_LOGGER.warn("stats: found {} reads without MateCigar attribute", MissingMateCigar);
        }
    }

    public void writeDuplicateStats(final ReduxConfig config)
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
            RD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }
}
