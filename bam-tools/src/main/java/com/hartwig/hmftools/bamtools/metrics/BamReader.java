package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.common.ReadGroup;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class BamReader
{
    private final MetricsConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final ChrBaseRegion mFilterRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final BaseCoverage mBaseCoverage;
    private final Map<String, ReadGroup> mReadGroupMap; // keyed by readId
    private final CombinedStats mCombinedStats;

    private int mTotalReads;
    private final PerformanceCounter mPerfCounter;
    private boolean mLogReadIds;

    public BamReader(
            final ChrBaseRegion region, final MetricsConfig config, final SamReader samReader, final BamSlicer bamSlicer,
            final CombinedStats combinedStats)
    {
        mConfig = config;
        mRegion = region;
        mCombinedStats = combinedStats;

        mSamReader = samReader;
        mBamSlicer = bamSlicer;

        mReadGroupMap = Maps.newHashMap();

        List<ChrBaseRegion> unmappableRegions = mConfig.UnmappableRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

        mBaseCoverage = new BaseCoverage(mConfig, mRegion.start(), mRegion.end(), unmappableRegions);

        mFilterRegion = null;

        mTotalReads = 0;
        mPerfCounter = new PerformanceCounter("Slice");

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        BT_LOGGER.debug("processing region({})", mRegion);

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

    private void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();

        ++mTotalReads;

        if(mTotalReads > 0 && (mTotalReads % 1_000_000) == 0)
            System.gc();

        if(mFilterRegion != null)
        {
            if(positionsOverlap(readStart, readStart + read.getReadBases().length, mFilterRegion.start(), mFilterRegion.end()))
                return;
        }

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(read.getReadName()))
                BT_LOGGER.debug("specific readId({}) unmapped({})", read.getReadName(), read.getReadUnmappedFlag());
        }

        // cache if the mate read overlaps
        ReadGroup readGroup = mReadGroupMap.get(read.getReadName());

        if(readGroup != null)
        {
            readGroup.addRead(read);

            if(readGroup.allReadsPresent())
            {
                processReadGroup(readGroup);
                mReadGroupMap.remove(read.getReadName());
            }

            return;
        }

        if(mateReadOverlaps(read))
        {
            readGroup = new ReadGroup(read);
            mReadGroupMap.put(readGroup.id(), readGroup);
            return;
        }

        // process this non-overlapping read immediately without caching
        mBaseCoverage.processRead(read, null);
    }

    private boolean mateReadOverlaps(final SAMRecord read)
    {
        if(mateUnmapped(read))
            return false;

        if(!read.getReferenceName().equals(read.getMateReferenceName()))
            return false;

        int readLength = read.getReadBases().length;
        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        int mateStart = read.getMateAlignmentStart();
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

        // check for overlaps in the reads and if found, form a list of the aligned coordinates
        boolean hasOverlaps = false;

        for(int i = 0; i < readGroup.reads().size() - 1; ++i)
        {
            SAMRecord read1 = readGroup.reads().get(i);

            for(int j = i + 1; j < readGroup.reads().size(); ++j)
            {
                SAMRecord read2 = readGroup.reads().get(j);

                if(positionsOverlap(read1.getAlignmentStart(), read1.getAlignmentEnd(), read2.getAlignmentStart(), read2.getAlignmentEnd()))
                {
                    hasOverlaps = true;
                    break;
                }
            }

            if(hasOverlaps)
                break;
        }

        if(!hasOverlaps)
        {
            readGroup.reads().forEach(x -> mBaseCoverage.processRead(x, null));
            return;
        }

        List<int[]> combinedAlignedBaseCoords = Lists.newArrayList();

        for(int i = 0; i < readGroup.reads().size(); ++i)
        {
            SAMRecord read = readGroup.reads().get(i);

            if(i == 0)
            {
                mBaseCoverage.processRead(read, null);
            }
            else
            {
                mBaseCoverage.processRead(read, combinedAlignedBaseCoords);
            }

            if(i < readGroup.reads().size() - 1)
                addAlignedCoords(read, combinedAlignedBaseCoords);
        }
    }

    private static void addAlignedCoords(final SAMRecord read, final List<int[]> alignedBaseCoords)
    {
        int position = read.getAlignmentStart();

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            switch(element.getOperator())
            {
                case S:
                case H:
                case I:
                    break;

                case M:
                    alignedBaseCoords.add(new int[] { position, position + element.getLength() - 1});
                    position += element.getLength();
                    break;

                case D:
                case N:
                    position += element.getLength();
                    break;

                default:
                    break;
            }
        }
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read)
    {
        processSamRecord(read);
    }

    @VisibleForTesting
    public BaseCoverage baseCoverage() { return mBaseCoverage; }

    @VisibleForTesting
    public Map<String,ReadGroup> readGroupMap() { return mReadGroupMap; }
}
