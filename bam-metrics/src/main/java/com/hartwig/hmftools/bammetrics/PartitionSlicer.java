package com.hartwig.hmftools.bammetrics;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bammetrics.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionSlicer
{
    private final BmConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final ChrBaseRegion mFilterRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final BaseCoverage mBaseCoverage;
    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId

    private final CombinedStats mCombinedStats;
    private int mTotalReads;
    private final PerformanceCounter mPerfCounter;

    private boolean mLogReadIds;

    public PartitionSlicer(
            final ChrBaseRegion region, final BmConfig config, final SamReader samReader, final BamSlicer bamSlicer,
            final CombinedStats combinedStats)
    {
        mConfig = config;
        mRegion = region;
        mCombinedStats = combinedStats;

        mSamReader = samReader;
        mBamSlicer = bamSlicer;

        mBaseCoverage = new BaseCoverage(mConfig, mRegion.start(), mRegion.end());

        mReadGroupMap = Maps.newHashMap();

        //ChrBaseRegion excludedRegion = getPolyGRegion(mConfig.RefGenVersion);
        // mFilterRegion = region.overlaps(excludedRegion) ? excludedRegion : null;
        mFilterRegion = null;

        mTotalReads = 0;
        mPerfCounter = new PerformanceCounter("Coverage");

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        BM_LOGGER.debug("processing region({})", mRegion);

        mPerfCounter.start(mConfig.PerfDebug ? mRegion.toString() : null);
        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        mPerfCounter.stop();

        // process overlapping groups
        for(ReadGroup readGroup : mReadGroupMap.values())
        {
            // determine overlapping bases and factor this into the coverage calcs
            processReadGroup(readGroup);
        }

        Metrics metrics = mBaseCoverage.createMetrics();

        mCombinedStats.addStats(metrics, mTotalReads, mPerfCounter);
    }

    private void processSamRecord(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();

        if(!mRegion.containsPosition(readStart))
            return;

        ++mTotalReads;

        if(mTotalReads > 0 && (mTotalReads % 1_000_000) == 0)
            System.gc();

        if(mFilterRegion != null)
        {
            if(positionsOverlap(readStart, readStart + record.getReadBases().length, mFilterRegion.start(), mFilterRegion.end()))
                return;
        }

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(record.getReadName()))
                BM_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
        }

        // cache if the mate read overlaps
        ReadGroup readGroup = mReadGroupMap.get(record.getReadName());

        if(readGroup != null)
        {
            readGroup.addRead(record);
            return;
        }

        if(mateReadOverlaps(record))
        {
            readGroup = new ReadGroup(record);
            mReadGroupMap.put(readGroup.id(), readGroup);
            return;
        }

        // process this non-overlapping read immediately without caching
        mBaseCoverage.processRead(record, null);
    }

    private boolean mateReadOverlaps(final SAMRecord record)
    {
        if(record.getMateUnmappedFlag())
            return false;

        if(!record.getReferenceName().equals(record.getMateReferenceName()))
            return false;

        int readLength = record.getReadBases().length;
        int readStart = record.getAlignmentStart();
        int readEnd = readStart + readLength - 1;

        int mateStart = record.getMateAlignmentStart();
        int mateEnd = mateStart + readLength - 1;

        return positionsOverlap(readStart, readEnd, mateStart, mateEnd);
    }

    private void processReadGroup(final ReadGroup readGroup)
    {
        if(readGroup.size() == 1)
        {
            mBaseCoverage.processRead(readGroup.reads().get(0), null);
            return;
        }

        List<int[]> readAlignedBaseCoords = Lists.newArrayList();

        for(int i = 0; i < readGroup.reads().size(); ++i)
        {
            SAMRecord record = readGroup.reads().get(i);

            int alignedBases = record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
            int[] alignedBaseCoords = new int[] { record.getAlignmentStart(), record.getAlignmentStart() + alignedBases - 1 };

            if(i == 0)
            {
                mBaseCoverage.processRead(record, null);
            }
            else
            {
                // check for any overlap with a previous read
                int[] overlappingBaseCoords = null;
                for(int j = 0; j < i; ++j)
                {
                    final int[] readBaseCoords = readAlignedBaseCoords.get(j);

                    if(positionsOverlap(alignedBaseCoords[SE_START], alignedBaseCoords[SE_END], readBaseCoords[SE_START], readBaseCoords[SE_END]))
                    {
                        overlappingBaseCoords = readBaseCoords;
                        break;
                    }
                }

                mBaseCoverage.processRead(record, overlappingBaseCoords);
            }

            readAlignedBaseCoords.add(alignedBaseCoords);
        }
    }
}
