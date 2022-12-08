package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.upperCoordsMatch;
import static com.hartwig.hmftools.bamtools.markdups.GroupCombiner.formChromosomePartition;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.BmConfig;
import com.hartwig.hmftools.bamtools.ReadGroup;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ChromosomeReader implements Consumer<PositionReadGroups>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final BaseRegion mCurrentPartition;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final GroupCombiner mGroupCombiner;
    private final RecordWriter mRecordWriter;
    private final ReadPositionArray mReadPositions;

    private final Map<String,SAMRecord> mSupplementaries;
    private final Map<String,ReadGroupInfo> mReadGroupInfos;

    private boolean mLogReadIds;

    public ChromosomeReader(
            final ChrBaseRegion region, final MarkDupsConfig config, final RecordWriter recordWriter,
            final GroupCombiner groupCombiner)
    {
        mConfig = config;
        mRegion = region;
        mGroupCombiner = groupCombiner;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
        // mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mRecordWriter = recordWriter;
        mReadPositions = new ReadPositionArray(region.Chromosome, config.BufferSize, this);

        mSupplementaries = Maps.newHashMap();
        mReadGroupInfos = Maps.newHashMap();
        mCurrentPartition = new BaseRegion(1, mConfig.PartitionSize);

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    @Override
    public Long call()
    {
        run();
        return (long)1; // return value not used
    }

    public void run()
    {
        BmConfig.BM_LOGGER.debug("processing chromosome({})", mRegion);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
    }

    private void processSamRecord(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();

        if(!mRegion.containsPosition(readStart))
            return;

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(record.getReadName()))
            {
                BmConfig.BM_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
            }
        }

        if(!mReadPositions.processRead(record))
        {
            // reasons for unhandled records:
            // - supplementaries not linked to a local group
            // - higher read with mate already processed (possibly in this partition)


            if(record.getSupplementaryAlignmentFlag())
            {
                mSupplementaries.put(record.getReadName(), record);
            }

        }

        // past the non-supplementary, lower reads to the array for initial duplicate evaluation, or other reads if are within the
        // current array bounds


        // otherwise
    }

    private void processSupplementary(final SAMRecord record)
    {
        SupplementaryReadData suppData = SupplementaryReadData.from(record);

        if(!suppData.Chromosome.equals(mRegion.Chromosome))
        {
            mGroupCombiner.handleSupplementary(record);
        }
        else if(suppData.Position > mCurrentPartition.end())
        {
            mSupplementaries.put(record.getReadName(), record);
        }
        else
        {
            // check dup status from earlier read groups
            ReadGroupInfo groupInfo = mReadGroupInfos.get(record.getReadName());

            if(groupInfo == null)
            {
                // a mistake has been made?
                mRecordWriter.writeRecord(record, false);
            }
            else
            {
                mRecordWriter.writeRecord(record, groupInfo.IsDuplicate);
            }
        }
    }

    public void accept(final PositionReadGroups positionReadGroups)
    {
        // identify duplicate groups, and cache any of these which aren't complete
        int i = 0;

        List<ReadGroup> readGroups = Lists.newArrayList(positionReadGroups.ReadGroups);
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();

        while(i < readGroups.size() - 1)
        {
            ReadGroup readGroup1 = readGroups.get(i);

            ChrBaseRegion lowerCoords1 = DuplicateGroup.lowerCoords(readGroup1);
            ChrBaseRegion upperCoords1 = DuplicateGroup.upperCoords(readGroup1);

            DuplicateGroup duplicateGroup = null;

            int j = i + 1;

            while(j < readGroups.size())
            {
                ReadGroup readGroup2 = readGroups.get(j);

                ChrBaseRegion lowerCoords2 = DuplicateGroup.lowerCoords(readGroup2);

                if(!lowerCoords1.matches(lowerCoords2))
                {
                    ++j;
                    continue;
                }

                ChrBaseRegion upperCoords2 = DuplicateGroup.upperCoords(readGroup2);

                if(!upperCoordsMatch(upperCoords1, upperCoords2))
                {
                    ++j;
                    continue;
                }

                // form a new group
                if(duplicateGroup == null)
                {
                    duplicateGroup = new DuplicateGroup(lowerCoords1, upperCoords1, readGroup1, readGroup2);
                    duplicateGroups.add(duplicateGroup);
                }
                else
                {
                    duplicateGroup.readGroups().add(readGroup2);
                }

                readGroups.remove(j);
            }

            if(duplicateGroup != null)
                readGroups.remove(i);
            else
                ++i;
        }
    }

    private class ReadGroupInfo
    {
        public final String ReadId;
        public final boolean IsDuplicate;
        public int MinPositionStart;
        public int MaxPositionStart;
        public final List<String> ChrPartitions;

        public ReadGroupInfo(final ReadGroup readGroup, boolean isDuplicate, final int partitionSize)
        {
            ReadId = readGroup.id();
            IsDuplicate = isDuplicate;

            List<String> chrPartitions = Lists.newArrayList();
            String chromosome = readGroup.reads().get(0).getContig();

            for(int i = 0; i < readGroup.reads().size(); ++i)
            {
                SAMRecord read = readGroup.reads().get(i);
                int readPosStart = read.getAlignmentStart();

                MinPositionStart = (i == 0) ? readPosStart : min(MinPositionStart, readPosStart);
                MaxPositionStart = max(MaxPositionStart, readPosStart);

                if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
                {
                    if(!read.getMateReferenceName().equals(chromosome))
                    {
                        chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                    }
                }

                SupplementaryReadData suppData = SupplementaryReadData.from(read);

                if(suppData != null && !suppData.Chromosome.equals(chromosome))
                    chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
            }

            ChrPartitions = !chrPartitions.isEmpty() ? chrPartitions : null;
        }

        /*
        public ReadGroupInfo(final String readId, final int minPositionStart, final boolean isDuplicate)
        {
            ReadId = readId;
            MinPositionStart = minPositionStart;
            IsDuplicate = isDuplicate;
            MaxPositionStart = 0;
        }
        */

        public String toString()
        {
            return format("positions(%d - %d) dup(%s) remotePartitions(%d) id(%s)",
                    MinPositionStart, MaxPositionStart, IsDuplicate, ChrPartitions != null ? ChrPartitions.size() : 0, ReadId);
        }
    }
}
