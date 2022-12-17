package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.calcBaseQualTotal;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.hasDuplicates;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroup.upperCoordsMatch;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.GroupCombiner.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.ReadGroupInfo.maxPositionStart;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.bamtools.ReadGroup;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
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

    private final List<SAMRecord> mPendingReads;
    private final List<DuplicateStatus> mPendingDuplicateState;

    private final boolean mLogReadIds;
    private int mTotalRecordCount;
    private final Set<SAMRecord> mReadsProcessed;
    private final DuplicateStats mStats;
    private final PerformanceCounter mPerfCounter;

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
        mTotalRecordCount = 0;

        mPendingReads = Lists.newArrayList();
        mPendingDuplicateState = Lists.newArrayList();

        mStats = new DuplicateStats();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mReadsProcessed = Sets.newHashSet();
        mPerfCounter = new PerformanceCounter("Slice");
    }

    public int totalRecordCount() { return mTotalRecordCount; }
    public Set<SAMRecord> readsProcessed() { return mReadsProcessed; }
    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public DuplicateStats duplicateStats() { return mStats; }

    @Override
    public Long call()
    {
        run();
        return (long)1;
    }

    public void run()
    {
        perfCounterStart();

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
            BM_LOGGER.info("processing chromosome({})", mRegion.Chromosome);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        }

        for(List<SAMRecord> supplementaries : mSupplementaries.values())
        {
            supplementaries.forEach(x -> mRecordWriter.writeRecord(x, DuplicateStatus.NONE));
        }

        onPartitionComplete(false);
        flushPendingReads();

        BM_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mTotalRecordCount);
        mSupplementaries.clear();
        mReadGroupInfos.clear();
        mReadsProcessed.clear();
    }

    private void onPartitionComplete(boolean setupNext)
    {
        mReadPositions.evictAll();
        mGroupCombiner.partitionComplete(mCurrentStrPartition);

        mPerfCounter.stop();

        BM_LOGGER.debug("partition({}:{}) complete, cached(supps={} readGroups={})",
                mRegion.Chromosome, mCurrentPartition, mSupplementaries.size(), mReadGroupInfos.size());

        if(setupNext)
        {
            mCurrentPartition.setStart(mCurrentPartition.end() + 1);
            mCurrentPartition.setEnd(mCurrentPartition.start() + mConfig.PartitionSize);
            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);

            // remove any group info with an expected start position now before the new partition
            Set<String> purgeReadIds = Sets.newHashSet();
            for(ReadGroupInfo groupInfo : mReadGroupInfos.values())
            {
                if(groupInfo.ExpectedRange.end() < mCurrentPartition.start())
                    purgeReadIds.add(groupInfo.ReadId);
            }

            purgeReadIds.forEach(x -> mReadGroupInfos.remove(x));

            perfCounterStart();
        }
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalRecordCount;

        if(BM_LOGGER.isTraceEnabled())
            mReadsProcessed.add(read);

        int readStart = read.getAlignmentStart();

        if(readStart > mCurrentPartition.end())
        {
            onPartitionComplete(true);
        }

        // if(!mCurrentPartition.containsPosition(readStart))
        //    return;

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(read.getReadName()))
            {
                BM_LOGGER.debug("specific readId({})", read.getReadName());
            }
        }

        if(!mReadPositions.processRead(read))
        {
            processNonPositionArrayRead(read);
        }
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
            writeRead(read, DuplicateStatus.NONE);
        }
        else
        {
            writeRead(read, groupInfo.Status);

            // remove this read group info if no further reads on this chromosome are expected
            int maxPositionStart = maxPositionStart(read);

            if(groupInfo.ExpectedRange.end() == maxPositionStart)
                mReadGroupInfos.remove(read.getReadName());
            else
                groupInfo.ExpectedRange.setEnd(maxPositionStart);
        }
    }

    private void processCandidateReadGroup(final ReadGroup readGroup, DuplicateStatus dupStatus)
    {
        ReadGroupInfo readGroupInfo = new ReadGroupInfo(readGroup, dupStatus, mCurrentPartition);

        writeReadGroup(readGroup, readGroupInfo.Status);

        if(readGroupInfo.IsComplete)
            return;

        if(readGroupInfo.ExpectedRange.start() < readGroupInfo.CurrentRange.start())
        {
            // check supplementaries
            List<SAMRecord> supplementaries = mSupplementaries.get(readGroup.id());

            if(supplementaries != null)
            {
                mSupplementaries.remove(readGroup.id());
                supplementaries.forEach(x -> writeRead(x, readGroupInfo.Status));
            }
        }

        if(readGroupInfo.ExpectedRange.end() > readGroupInfo.CurrentRange.end())
        {
            // mate reads in a later partition expected
            mReadGroupInfos.put(readGroup.id(), readGroupInfo);
        }

        if(readGroupInfo.ChrPartitions != null)
        {
            mGroupCombiner.processReadGroup(readGroup.id(), readGroupInfo.ChrPartitions, readGroupInfo.Status);
        }
    }

    public void accept(final PositionReadGroups positionReadGroups)
    {
        // identify duplicate groups, and cache any of these which aren't complete
        int i = 0;

        List<ReadGroup> readGroups = positionReadGroups.ReadGroups.values().stream().collect(Collectors.toList());
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
            ReadGroup primaryGroup = duplicateGroup.findPrimaryGroup(true);

            if(mConfig.RunChecks)
            {
                // log discrepancies
                ReadGroup calcPrimaryGroup = duplicateGroup.findPrimaryGroup(false);

                boolean logDiscrepancy = false;

                if(primaryGroup != calcPrimaryGroup && calcBaseQualTotal(primaryGroup) != calcBaseQualTotal(calcPrimaryGroup))
                    logDiscrepancy = true;
                else if(duplicateGroup.readGroups().stream().anyMatch(x -> hasDuplicates(x) == (x == primaryGroup)))
                    logDiscrepancy = true;

                if(logDiscrepancy)
                {
                    for(ReadGroup readGroup : duplicateGroup.readGroups())
                    {
                        BM_LOGGER.trace("readGroup({}) hasDups({}) isPrimary({}) baseQualTotal({})",
                                readGroup.toString(), hasDuplicates(readGroup), readGroup == primaryGroup, calcBaseQualTotal(readGroup));
                    }
                }
            }

            for(ReadGroup readGroup : duplicateGroup.readGroups())
            {
                DuplicateStatus dupStatus = readGroup == primaryGroup ? PRIMARY : DUPLICATE;
                processCandidateReadGroup(readGroup, dupStatus);
            }
        }

        readGroups.forEach(x -> processCandidateReadGroup(x, DuplicateStatus.NONE));
    }

    private void writeReadGroup(final ReadGroup readGroup, DuplicateStatus duplicateStatus)
    {
        readGroup.reads().forEach(x -> writeRead(x, duplicateStatus));
    }

    private void writeRead(final SAMRecord read, DuplicateStatus duplicateStatus)
    {
        if(mConfig.WriteCacheSize == 0)
        {
            mRecordWriter.writeRecord(read, duplicateStatus);
            return;
        }

        if(mPendingReads.size() >= mConfig.WriteCacheSize)
        {
            flushPendingReads();
            mRecordWriter.writeRecord(read, duplicateStatus);
        }
        else
        {
            mPendingReads.add(read);
            mPendingDuplicateState.add(duplicateStatus);
        }
    }

    private void flushPendingReads()
    {
        for(int i = 0; i < mPendingReads.size(); ++i)
        {
            mRecordWriter.writeRecord(mPendingReads.get(i), mPendingDuplicateState.get(i));
        }

        mPendingReads.clear();
        mPendingDuplicateState.clear();
    }

    private void perfCounterStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPerfCounter.start();
    }
}
