package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.consensus.UmiUtils.calcUmiIdDiff;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.MarkDupsConfig;
import com.hartwig.hmftools.markdups.consensus.UmiConfig;
import com.hartwig.hmftools.markdups.consensus.DuplicateGroup;

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

    public final Map<Integer,Integer> DuplicateFrequencies;
    public final List<UmiGroupCounts> UmiGroupFrequencies;
    public final List<PositionFragmentsData> PositionFragments;

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
        UmiGroupFrequencies = Lists.newArrayList();
        PositionFragments = Lists.newArrayList();
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

        for(UmiGroupCounts otherUgCounts : other.UmiGroupFrequencies)
        {
            UmiGroupCounts matchedUgCounts = UmiGroupFrequencies.stream()
                    .filter(x -> x.DuplicateGroupCount == otherUgCounts.DuplicateGroupCount
                    && x.ReadCount == otherUgCounts.ReadCount)
                    .findFirst().orElse(null);

            if(matchedUgCounts != null)
            {
                matchedUgCounts.GroupCount += otherUgCounts.GroupCount;

                for(int i = 0; i < otherUgCounts.EditDistanceFrequency.length; ++i)
                {
                    matchedUgCounts.EditDistanceFrequency[i] += otherUgCounts.EditDistanceFrequency[i];
                }
            }
            else
            {
                UmiGroupFrequencies.add(otherUgCounts);
            }
        }

        for(PositionFragmentsData otherPosFragData : other.PositionFragments)
        {
            PositionFragmentsData posFragData = getOrCreatePositionFragmentData(otherPosFragData.PosGroupCount, otherPosFragData.UniqueFragmentCount);

            if(posFragData.Frequency == 0)
            {
                posFragData.UmiGroupDetails = otherPosFragData.UmiGroupDetails;
                posFragData.MaxPosUmiCount = otherPosFragData.MaxPosUmiCount;
                posFragData.MaxUmiReadsCount = otherPosFragData.MaxUmiReadsCount;
            }

            posFragData.Frequency += otherPosFragData.Frequency;

        }
    }

    private static int roundFrequency(int frequency)
    {
        if(frequency <= 10)
            return frequency;
        else if(frequency <= 100)
            return round(frequency/10) * 10;
        else if(frequency <= 1000)
            return round(frequency/100) * 100;
        else
            return round(frequency/1000) * 1000;
    }

    public void addFrequency(int frequency)
    {
        int rounded = roundFrequency(frequency);
        Integer count = DuplicateFrequencies.get(rounded);
        DuplicateFrequencies.put(rounded, count == null ? 1 : count + 1);
    }

    private class PositionFragmentsData
    {
        public final int PosGroupCount;
        public final int UniqueFragmentCount;
        public int Frequency;
        public int MaxPosUmiCount;
        public int MaxUmiReadsCount; // max reads collapsed into a UMI group
        public String UmiGroupDetails;

        public PositionFragmentsData(final int posGroupCount, final int uniqueFragmentCount)
        {
            PosGroupCount = posGroupCount;
            UniqueFragmentCount = uniqueFragmentCount;
            Frequency = 0;
            MaxPosUmiCount = 0;
            MaxUmiReadsCount = 0;
            UmiGroupDetails = "";
        }

        public String toString()
        {
            return format("pos(%d) unique(%d) count(%d) max(umi=%d reads=%d)",
                    PosGroupCount, UniqueFragmentCount, Frequency, MaxPosUmiCount, MaxUmiReadsCount);
        }
    }

    private PositionFragmentsData getOrCreatePositionFragmentData(int posGroupCount, int duplicatePosCount)
    {
        PositionFragmentsData matchedPosFragments = PositionFragments.stream()
                .filter(x -> x.PosGroupCount == posGroupCount && x.UniqueFragmentCount == duplicatePosCount).findFirst().orElse(null);

        if(matchedPosFragments == null)
        {
            matchedPosFragments = new PositionFragmentsData(posGroupCount, duplicatePosCount);
            PositionFragments.add(matchedPosFragments);
        }

        return matchedPosFragments;
    }

    public void recordFragmentPositions(int posGroupCount, int duplicatePosCount, int maxDuplicatePosCount, final DuplicateGroup duplicateGroup)
    {
        PositionFragmentsData posFragData = getOrCreatePositionFragmentData(posGroupCount, duplicatePosCount);
        ++posFragData.Frequency;
        posFragData.MaxPosUmiCount = max(posFragData.MaxPosUmiCount, maxDuplicatePosCount);

        if(duplicateGroup != null && duplicateGroup.fragmentCount() > posFragData.MaxUmiReadsCount)
        {
            posFragData.MaxUmiReadsCount = duplicateGroup.fragmentCount();
            posFragData.UmiGroupDetails = format("%s %s", duplicateGroup.coordinatesKey(), duplicateGroup.getReadIds().get(0));
        }
    }

    private static final int MAX_EDIT_DISTANCE = 10;

    private class UmiGroupCounts
    {
        public final int DuplicateGroupCount;
        public final int ReadCount;
        public int GroupCount;

        public final int[] EditDistanceFrequency;

        public UmiGroupCounts(final int duplicateGroupCount, final int readCount)
        {
            DuplicateGroupCount = duplicateGroupCount;
            ReadCount = readCount;
            GroupCount = 0;
            EditDistanceFrequency = new int[MAX_EDIT_DISTANCE+1];
        }
    }

    public void recordUmiBaseDiffs(final UmiConfig umiConfig, final List<DuplicateGroup> umiGroups)
    {
        // evaluate 1 or 2 UMI groups, including those with a single fragment which may have been under-clustered
        if(umiGroups.size() == 1)
        {
            recordUmiGroupStats(umiConfig, umiGroups.get(0));
        }
        else if(umiGroups.size() == 2)
        {
            recordUmiGroupStats(umiConfig, umiGroups.get(0), umiGroups.get(1));
        }
    }

    private void recordUmiGroupStats(final UmiConfig umiConfig, final DuplicateGroup umiGroup)
    {
        UmiGroupCounts umiGroupStats = getOrCreateUmiGroupCounts(1, umiGroup.fragmentCount());
        ++umiGroupStats.GroupCount;

        for(String readId : umiGroup.getReadIds())
        {
            int diff = calcUmiIdDiff(umiConfig.extractUmiId(readId), umiGroup.id());

            if(diff <= MAX_EDIT_DISTANCE)
                ++umiGroupStats.EditDistanceFrequency[diff];
        }
    }

    private void recordUmiGroupStats(final UmiConfig umiConfig, final DuplicateGroup group1, final DuplicateGroup group2)
    {
        UmiGroupCounts umiGroupStats = getOrCreateUmiGroupCounts(2, group1.fragmentCount() + group2.fragmentCount());
        ++umiGroupStats.GroupCount;

        for(int groupIndex = 0; groupIndex <= 1; ++groupIndex)
        {
            DuplicateGroup testGroup = groupIndex == 0 ? group1 : group2;
            DuplicateGroup readsGroup = groupIndex == 0 ? group2 : group1;

            for(String readId : readsGroup.getReadIds())
            {
                int diff = calcUmiIdDiff(umiConfig.extractUmiId(readId), testGroup.id());

                if(diff <= MAX_EDIT_DISTANCE)
                    ++umiGroupStats.EditDistanceFrequency[diff];
            }
        }
    }

    private UmiGroupCounts getOrCreateUmiGroupCounts(int groupCount, int readCount)
    {
        int roundedReadCount = roundFrequency(readCount);

        UmiGroupCounts matchedStats = UmiGroupFrequencies.stream()
                .filter(x -> x.DuplicateGroupCount == groupCount && x.ReadCount == roundedReadCount)
                .findFirst().orElse(null);

        if(matchedStats == null)
        {
            matchedStats = new UmiGroupCounts(groupCount, roundedReadCount);
            UmiGroupFrequencies.add(matchedStats);
        }

        return matchedStats;
    }

    public void addDuplicateGroup(final int fragmentCount)
    {
        addFrequency(fragmentCount);
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
                MD_LOGGER.debug("duplicate frequency({}={})", frequency, DuplicateFrequencies.get(frequency));
            }
        }
    }

    public void writeDuplicateStats(final MarkDupsConfig config)
    {
        try
        {
            String filename = config.formFilename("duplicate_stats");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("DuplicateReadCount,Frequency");
            writer.newLine();

            List<Integer> frequencies = DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            for(Integer frequency : frequencies)
            {
                writer.write(format("%d,%d", frequency, DuplicateFrequencies.get(frequency)));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }

    public void writeUmiStats(final MarkDupsConfig config)
    {
        try
        {
            String filename = config.formFilename("umi_stats");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("DuplicateGroupCount,ReadCount,GroupCount");

            for(int i = 0; i <= MAX_EDIT_DISTANCE; ++i)
            {
                writer.write(format(",ED_%d", i));
            }

            writer.newLine();

            for(UmiGroupCounts umiGroupCounts : UmiGroupFrequencies)
            {
                writer.write(format("%d,%d,%d", umiGroupCounts.DuplicateGroupCount, umiGroupCounts.ReadCount, umiGroupCounts.GroupCount));

                for(int i = 0; i <= MAX_EDIT_DISTANCE; ++i)
                {
                    writer.write(format(",%d", umiGroupCounts.EditDistanceFrequency[i]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }

    public void writePositionFragmentsData(final MarkDupsConfig config)
    {
        try
        {
            String filename = config.formFilename("pos_frags_stats");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("UniqueFragCoords,UniqueFragments,Count,MaxUmis,MaxUmiDuplicates,MaxUmiDetails");
            writer.newLine();

            for(PositionFragmentsData posFragData : PositionFragments)
            {
                writer.write(format("%d,%d,%d,%d,%d,%s",
                        posFragData.PosGroupCount, posFragData.UniqueFragmentCount, posFragData.Frequency,
                        posFragData.MaxPosUmiCount, posFragData.MaxUmiReadsCount, posFragData.UmiGroupDetails));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write position fragments stats: {}", e.toString());
        }
    }
}
