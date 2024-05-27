package com.hartwig.hmftools.redux.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.umi.UmiUtils.calcUmiIdDiff;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.DuplicateGroup;

public class UmiStatistics
{
    public long UmiGroups;
    public final List<UmiGroupCounts> UmiGroupFrequencies;
    public final List<PositionFragmentCounts> PositionFragments;
    public int[][] UmiPositionBaseFrequencies;

    public static final int UMI_BASE_COUNT = Nucleotides.DNA_BASES.length + 1; // will record any Ns

    public UmiStatistics()
    {
        UmiGroups = 0;
        UmiGroupFrequencies = Lists.newArrayList();
        PositionFragments = Lists.newArrayList();
        UmiPositionBaseFrequencies = null;
    }

    public void merge(final UmiStatistics other)
    {
        UmiGroups += other.UmiGroups;

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

        for(PositionFragmentCounts otherPosFragData : other.PositionFragments)
        {
            PositionFragmentCounts posFragData = getOrCreatePositionFragmentData(otherPosFragData.UniqueCoordCount, otherPosFragData.UniqueFragmentCount);

            if(posFragData.Frequency == 0)
            {
                posFragData.UmiGroupDetails = otherPosFragData.UmiGroupDetails;
                posFragData.MaxCoordUmiCount = otherPosFragData.MaxCoordUmiCount;
                posFragData.MaxUmiReadsCount = otherPosFragData.MaxUmiReadsCount;
            }

            posFragData.Frequency += otherPosFragData.Frequency;

        }

        if(other.UmiPositionBaseFrequencies != null)
        {
            if(UmiPositionBaseFrequencies == null)
            {
                UmiPositionBaseFrequencies = new int[other.UmiPositionBaseFrequencies.length][UMI_BASE_COUNT];
            }

            for(int p = 0; p < UmiPositionBaseFrequencies.length; ++p)
            {
                for(int b = 0; b < UMI_BASE_COUNT; ++b)
                {
                    UmiPositionBaseFrequencies[p][b] += other.UmiPositionBaseFrequencies[p][b];
                }
            }
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

    private PositionFragmentCounts getOrCreatePositionFragmentData(int uniqueCoordCount, int uniqueFragmentCount)
    {
        PositionFragmentCounts matchedPosFragments = PositionFragments.stream()
                .filter(x -> x.UniqueCoordCount == uniqueCoordCount && x.UniqueFragmentCount == uniqueFragmentCount).findFirst().orElse(null);

        if(matchedPosFragments == null)
        {
            matchedPosFragments = new PositionFragmentCounts(uniqueCoordCount, uniqueFragmentCount);
            PositionFragments.add(matchedPosFragments);
        }

        return matchedPosFragments;
    }

    public void recordFragmentPositions(int uniqueCoordCount, int uniqueFragmentCount, int maxCoordUmiCount, final DuplicateGroup duplicateGroup)
    {
        PositionFragmentCounts posFragData = getOrCreatePositionFragmentData(uniqueCoordCount, uniqueFragmentCount);
        ++posFragData.Frequency;
        posFragData.MaxCoordUmiCount = max(posFragData.MaxCoordUmiCount, maxCoordUmiCount);

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

    public void recordUmiBaseStats(final UmiConfig umiConfig, final List<DuplicateGroup> umiGroups)
    {
        umiGroups.forEach(x -> recordUmiBaseFrequencies(x));

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

    private void recordUmiBaseFrequencies(final DuplicateGroup umiGroup)
    {
        String umiId = umiGroup.umiId();

        if(UmiPositionBaseFrequencies == null)
            UmiPositionBaseFrequencies = new int[umiId.length()][UMI_BASE_COUNT];

        for(int p = 0; p < min(umiId.length(), UmiPositionBaseFrequencies.length); ++p)
        {
            int baseIndex = getBaseIndex(umiGroup.umiId().charAt(p));

            if(baseIndex >= 0)
                ++UmiPositionBaseFrequencies[p][baseIndex];
        }
    }

    private static int getBaseIndex(final char base)
    {
        switch(base)
        {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            case 'N': return 4;
        }

        return -1;
    }

    private void recordUmiGroupStats(final UmiConfig umiConfig, final DuplicateGroup umiGroup)
    {
        UmiGroupCounts umiGroupStats = getOrCreateUmiGroupCounts(1, umiGroup.fragmentCount());
        ++umiGroupStats.GroupCount;

        for(String readId : umiGroup.getReadIds())
        {
            int diff = calcUmiIdDiff(umiConfig.extractUmiId(readId), umiGroup.umiId());

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
                int diff = calcUmiIdDiff(umiConfig.extractUmiId(readId), testGroup.umiId());

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

    public void writeUmiBaseDiffStats(final ReduxConfig config)
    {
        try
        {
            String filename = config.formFilename("umi_edit_distance");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("UmiCountWithCoords\tDuplicateReadCount\tFrequency");

            for(int i = 0; i <= MAX_EDIT_DISTANCE; ++i)
            {
                writer.write(format("\tED_%d", i));
            }

            writer.newLine();

            for(UmiGroupCounts umiGroupCounts : UmiGroupFrequencies)
            {
                writer.write(format("%d\t%d\t%d", umiGroupCounts.DuplicateGroupCount, umiGroupCounts.ReadCount, umiGroupCounts.GroupCount));

                for(int i = 0; i <= MAX_EDIT_DISTANCE; ++i)
                {
                    writer.write(format("\t%d", umiGroupCounts.EditDistanceFrequency[i]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }

    public void writeUmiBaseFrequencyStats(final ReduxConfig config)
    {
        try
        {
            String filename = config.formFilename("umi_nucleotide_freq");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("UmiPosition\tACount\tCCount\tGCount\tTCount\tNCount");
            writer.newLine();

            for(int p = 0; p < UmiPositionBaseFrequencies.length; ++p)
            {
                writer.write(format("%d", p+1));

                for(int b = 0; b < UMI_BASE_COUNT; ++b)
                {
                    writer.write(format("\t%d", UmiPositionBaseFrequencies[p][b]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write UMI stats: {}", e.toString());
        }
    }
    public void writePositionFragmentsData(final ReduxConfig config)
    {
        try
        {
            String filename = config.formFilename("umi_coord_freq");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("UniqueCoordsWithStartPos\tUniquePrimariesWithStartPos\tFrequency\tMaxCoordUmiCount\tMaxUmiReads\tMaxUmiDetails");
            writer.newLine();

            for(PositionFragmentCounts posFragData : PositionFragments)
            {
                writer.write(format("%d\t%d\t%d\t%d\t%d\t%s",
                        posFragData.UniqueCoordCount, posFragData.UniqueFragmentCount, posFragData.Frequency,
                        posFragData.MaxCoordUmiCount, posFragData.MaxUmiReadsCount, posFragData.UmiGroupDetails));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write position fragments stats: {}", e.toString());
        }
    }
}
