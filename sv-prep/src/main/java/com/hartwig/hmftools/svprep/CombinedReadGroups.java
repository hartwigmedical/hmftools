package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

public class CombinedReadGroups
{
    private final Map<String,Map<String,ReadGroup>> mIncompleteReadGroups; // keyed by chromosome then readId
    private final PerformanceCounter mPerfCounter;

    public CombinedReadGroups()
    {
        mIncompleteReadGroups = Maps.newHashMap();
        mPerfCounter = new PerformanceCounter("ReadMerge");
    }

    public void logPerfStats()
    {
        mPerfCounter.logStats();
    }

    public synchronized List<ReadGroup> addIncompleteReadGroup(
            final String chrPartition, final Map<String,Map<String,ReadGroup>> chrIncompleteGroups)
    {
        List<ReadGroup> completeGroups = Lists.newArrayList();

        if(chrIncompleteGroups.isEmpty())
            return completeGroups;

        mPerfCounter.start();

        int initTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();
        int initChrIncomplete = chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum();

        // incomplete groups are looked up by the new chromosome partition to check for matches, and any which aren't found and
        // stored against the remote chromosome partition

        for(Map.Entry<String,Map<String, ReadGroup>> entry : chrIncompleteGroups.entrySet())
        {
            String otherChrPartition = entry.getKey();
            Map<String, ReadGroup> newIncompleteGroups = entry.getValue();

            Map<String,ReadGroup> existingGroups = mIncompleteReadGroups.get(otherChrPartition);

            if(existingGroups == null)
            {
                Map<String,ReadGroup> chrPartitionGroups = mIncompleteReadGroups.get(chrPartition);
                if(chrPartitionGroups == null)
                {
                    chrPartitionGroups = Maps.newHashMap();
                    mIncompleteReadGroups.put(chrPartition, chrPartitionGroups);
                }

                chrPartitionGroups.putAll(newIncompleteGroups);

                SV_LOGGER.trace("added chromosome partitions pair({} & {}) newGroups({})",
                        chrPartition, otherChrPartition, newIncompleteGroups.size());
            }
            else
            {
                mergeReadMaps(existingGroups, completeGroups, newIncompleteGroups);

                SV_LOGGER.trace("combined chromosome partitions pair({} & {}) existing({}) new({}) complete({})",
                        chrPartition, otherChrPartition, existingGroups.size(), newIncompleteGroups.size(), completeGroups.size());
            }
        }

        int newTotalIncomplete = mIncompleteReadGroups.values().stream().mapToInt(x -> x.size()).sum();

        SV_LOGGER.info("chromosomePartition({}) complete({}) partials chrPartition({}) total({} -> {})",
                chrPartition, completeGroups.size(), initChrIncomplete, initTotalIncomplete, newTotalIncomplete);

        mPerfCounter.stop();

        return completeGroups;
    }

    public void writeRemainingReadGroups(final ResultsWriter writer, final Set<WriteType> writeTypes)
    {
        if(!writeTypes.contains(BAM) && !writeTypes.contains(READS))
            return;

        Map<String,ReadGroup> readGroups = Maps.newHashMap();

        for(Map<String,ReadGroup> readGroupMaps : mIncompleteReadGroups.values())
        {
            for(ReadGroup readGroup : readGroupMaps.values())
            {
                ReadGroup existingGroup = readGroups.get(readGroup.id());

                if(existingGroup != null)
                {
                    existingGroup.merge(readGroup);
                    existingGroup.setGroupStatus(null);
                }
                else
                {
                    readGroups.put(readGroup.id(), readGroup);
                }
            }
        }

        if(writeTypes.contains(BAM))
            readGroups.values().forEach(x -> writer.writeBamRecords(x));

        if(writeTypes.contains(READS))
            writer.writeReadData(readGroups.values().stream().collect(Collectors.toList()));

        if(!readGroups.isEmpty())
        {
            SV_LOGGER.info("remaining partial groups({})", readGroups.size());
        }

        mPerfCounter.logStats();
    }

    private static final String CHR_PARTITION_DELIM = "_";

    public static String externalReadChrPartition(final ChrBaseRegion region, int partitionSize, final List<ReadRecord> reads)
    {
        for(ReadRecord read : reads)
        {
            if(!region.containsPosition(read.MateChromosome, read.MatePosStart))
            {
                return formChromosomePartition(read.MateChromosome, read.MatePosStart, partitionSize);
            }
            else if(read.hasSuppAlignment())
            {
                return formChromosomePartition(
                        read.supplementaryAlignment().Chromosome, read.supplementaryAlignment().Position, partitionSize);
            }
        }

        return null;
    }

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    private void mergeReadMaps(
            final Map<String,ReadGroup> partialGroups, final List<ReadGroup> completeGroups, final Map<String,ReadGroup> sourceMap)
    {
        // 1. copies complete groups from the source map into complete groups map
        // 2. checks for a partial match by combining partials and source, and if found removes from partials
        // 3. new partial groups from the source map are copied into the partials map
        // note: source map is logically const
        for(Map.Entry<String, ReadGroup> entry : sourceMap.entrySet())
        {
            final ReadGroup srcReadGroup = entry.getValue();

            if(srcReadGroup.isComplete())
            {
                completeGroups.add(srcReadGroup);
            }
            else
            {
                // look for an existing incomplete group to add these reads to
                final String readId = entry.getKey();
                ReadGroup existingReadGroup = partialGroups.get(readId);

                if(existingReadGroup == null)
                {
                    partialGroups.put(readId, srcReadGroup);
                }
                else
                {
                    existingReadGroup.merge(srcReadGroup);

                    if(existingReadGroup.isComplete())
                    {
                        partialGroups.remove(readId);
                        completeGroups.add(existingReadGroup);
                    }
                }
            }
        }
    }
}
