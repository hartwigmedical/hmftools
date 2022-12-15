package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.calcBaseQualTotal;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.hasDuplicates;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.upperCoordsMatch;
import static com.hartwig.hmftools.bamtools.markdups.GroupCombiner.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.ReadGroupInfo.maxPositionStart;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
    private String mCurrentStrPartition;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final GroupCombiner mGroupCombiner;
    private final RecordWriter mRecordWriter;
    private final ReadPositionArray mReadPositions;

    private final Map<String,List<SAMRecord>> mSupplementaries;
    private final Map<String,ReadGroupInfo> mReadGroupInfos;

    private final boolean mLogReadIds;

    public ChromosomeReader(
            final ChrBaseRegion region, final MarkDupsConfig config, final RecordWriter recordWriter,
            final GroupCombiner groupCombiner)
    {
        mConfig = config;
        mRegion = region;
        mGroupCombiner = groupCombiner;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        // mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mRecordWriter = recordWriter;
        mReadPositions = new ReadPositionArray(region.Chromosome, config.BufferSize, this);

        mSupplementaries = Maps.newHashMap();
        mReadGroupInfos = Maps.newHashMap();

        if(!mConfig.SpecificRegions.isEmpty())
        {
            ChrBaseRegion firstRegion = mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(mRegion.Chromosome)).findFirst().orElse(mRegion);
            int partitionStart = (firstRegion.start() / mConfig.PartitionSize) * mConfig.PartitionSize;
            mCurrentPartition = new BaseRegion(partitionStart, partitionStart + mConfig.PartitionSize - 1);
        }
        else
        {
            mCurrentPartition = new BaseRegion(1, mConfig.PartitionSize);
        }

        mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    @Override
    public Long call()
    {
        run();
        return (long)1;
    }

    public void run()
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            for(ChrBaseRegion region : mConfig.SpecificRegions)
            {
                if(!region.Chromosome.equals(mRegion.Chromosome))
                    continue;

                BM_LOGGER.debug("processing specific region({})", region);
                mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processSamRecord);
            }
        }
        else
        {
            BM_LOGGER.debug("processing chromosome({})", mRegion);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        }

        mReadPositions.evictAll();
    }

    private void processSamRecord(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();

        if(readStart > mCurrentPartition.end())
        {
            setNextPartition();
        }

        // if(!mCurrentPartition.containsPosition(readStart))
        //    return;

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(record.getReadName()))
            {
                BM_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
            }
        }


        if(!mReadPositions.processRead(record))
        {
            processNonPositionArrayRead(record);
        }
    }

    private void setNextPartition()
    {
        mGroupCombiner.partitionComplete(mCurrentStrPartition);
        mCurrentPartition.setStart(mCurrentPartition.end() + 1);
        mCurrentPartition.setEnd(mCurrentPartition.start() + mConfig.PartitionSize);
        mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
        mReadPositions.evictAll();
    }

    private void processNonPositionArrayRead(final SAMRecord read)
    {
        // reasons for unhandled records:
        // - supplementaries not linked to a local group
        // - higher read with mate already processed (possibly in this partition)

        // check if this record's group has been handled already in this partition
        if(read.getSupplementaryAlignmentFlag())
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read);

            if(!suppData.Chromosome.equals(mRegion.Chromosome) || !mCurrentPartition.containsPosition(suppData.Position))
            {
                mGroupCombiner.handleSupplementary(read, mCurrentStrPartition);
                return;
            }
            else if(suppData.Position > read.getAlignmentStart())
            {
                List<SAMRecord> supplementaries = mSupplementaries.get(read.getReadName());

                if(supplementaries == null)
                {
                    supplementaries = Lists.newArrayList();
                    mSupplementaries.put(read.getReadName(), supplementaries);
                }

                supplementaries.add(read);
                return;
            }
        }

        // check for a previously processed group
        ReadGroupInfo groupInfo = mReadGroupInfos.get(read.getReadName());

        if(groupInfo == null)
        {
            // a mistake has been made?
            mRecordWriter.writeRecord(read, false);
        }
        else
        {
            mRecordWriter.writeRecord(read, groupInfo.IsDuplicate);

            // remove this read group info if no further reads on this chromosome are expected
            int maxPositionStart = maxPositionStart(read);

            if(groupInfo.ExpectedRange.end() == maxPositionStart)
                mReadGroupInfos.remove(read.getReadName());
            else
                groupInfo.ExpectedRange.setEnd(maxPositionStart);
        }
    }

    private void processCandidateReadGroup(final ReadGroup readGroup, boolean isDuplicate)
    {
        ReadGroupInfo readGroupInfo = new ReadGroupInfo(readGroup, isDuplicate, mCurrentPartition);

        mRecordWriter.writeReadGroup(readGroup, readGroupInfo.IsDuplicate);

        if(readGroupInfo.IsComplete)
            return;

        if(readGroupInfo.ExpectedRange.start() < readGroupInfo.CurrentRange.start())
        {
            // check supplementaries
            List<SAMRecord> supplementaries = mSupplementaries.get(readGroup.id());

            if(supplementaries != null)
            {
                mSupplementaries.remove(readGroup.id());
                supplementaries.forEach(x -> mRecordWriter.writeRecord(x, readGroupInfo.IsDuplicate));
            }
        }

        if(readGroupInfo.ExpectedRange.end() > readGroupInfo.CurrentRange.end())
        {
            // mate reads in a later partition expected
            mReadGroupInfos.put(readGroup.id(), readGroupInfo);
        }

        if(readGroupInfo.ChrPartitions != null)
        {
            mGroupCombiner.processReadGroup(readGroup.id(), readGroupInfo.ChrPartitions, readGroupInfo.IsDuplicate);
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

        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            // primary non-duplicate is the fragment with the highest base qual
            ReadGroup primaryGroup = duplicateGroup.findPrimaryGroup();

            if(duplicateGroup.readGroups().stream().anyMatch(x -> hasDuplicates(x) == (x == primaryGroup)))
            {
                for(ReadGroup readGroup : duplicateGroup.readGroups())
                {
                    BM_LOGGER.debug("readGroup({}) hasDups({}) isPrimary({}) baseQualTotal({})",
                            readGroup.toString(), hasDuplicates(readGroup), readGroup == primaryGroup, calcBaseQualTotal(readGroup));
                }
            }

            for(ReadGroup readGroup : duplicateGroup.readGroups())
            {
                processCandidateReadGroup(readGroup, readGroup != primaryGroup);
            }
        }

        readGroups.forEach(x -> processCandidateReadGroup(x, false));
    }
}
