package com.hartwig.hmftools.bamtools;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.metrics.BaseCoverage;
import com.hartwig.hmftools.bamtools.metrics.CombinedStats;
import com.hartwig.hmftools.bamtools.metrics.Metrics;
import com.hartwig.hmftools.bamtools.slice.SliceWriter;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
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

    // coverage
    private final boolean mRunMetrics;
    private final BaseCoverage mBaseCoverage;
    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final CombinedStats mCombinedStats;

    // slicing
    private final boolean mRunSlice;
    private final List<BaseRegion> mSliceRegions;
    private final SliceWriter mSliceWriter;

    private int mTotalReads;
    private final PerformanceCounter mPerfCounter;
    private boolean mLogReadIds;

    public PartitionSlicer(
            final ChrBaseRegion region, final BmConfig config, final SamReader samReader, final BamSlicer bamSlicer,
            final CombinedStats combinedStats, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mRegion = region;
        mCombinedStats = combinedStats;
        mSliceWriter = sliceWriter;

        mSamReader = samReader;
        mBamSlicer = bamSlicer;

        mRunMetrics = config.runMetrics();
        mRunSlice = config.runSlicing();

        mReadGroupMap = Maps.newHashMap();

        mBaseCoverage = mRunMetrics ? new BaseCoverage(mConfig, mRegion.start(), mRegion.end()) : null;

        mSliceRegions = mConfig.SpecificRegions.stream()
                .filter(x -> mRegion.chromosome().equals(x.chromosome()))
                .filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                .map(x -> new BaseRegion(x.start(), x.end()))
                .collect(Collectors.toList());

        mFilterRegion = null;

        mTotalReads = 0;
        mPerfCounter = new PerformanceCounter("Slice");

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        BmConfig.BM_LOGGER.debug("processing region({})", mRegion);

        mPerfCounter.start(mConfig.PerfDebug ? mRegion.toString() : null);
        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        mPerfCounter.stop();

        if(mConfig.runMetrics())
        {
            // process overlapping groups
            for(ReadGroup readGroup : mReadGroupMap.values())
            {
                // determine overlapping bases and factor this into the coverage calcs
                processReadGroup(readGroup);
            }

            Metrics metrics = mBaseCoverage.createMetrics();
            mCombinedStats.addStats(metrics, mTotalReads, mPerfCounter);
        }
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
                BmConfig.BM_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
        }

        if(mRunSlice)
            handleSliceRecord(record);

        if(mRunMetrics)
            handleMetricsRecord(record);
    }

    private void handleSliceRecord(final SAMRecord record)
    {
        // any read which overlaps a slicing region or its mate does, or its supplementary does

        boolean overlapsSliceRegion = false;
        int readLength = record.getBaseQualities().length - 1;
        boolean hasRemoteSupplementary = false;

        if(overlapsSliceRegion(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd()))
        {
            overlapsSliceRegion = true;
        }
        else if(overlapsSliceRegion(
                record.getMateReferenceName(), record.getMateAlignmentStart(), record.getMateAlignmentStart() + readLength))
        {
            overlapsSliceRegion = true;
        }
        else if(record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            if(suppData != null)
            {
                if(overlapsSliceRegion(suppData.Chromosome, suppData.Position, suppData.Position + readLength))
                    overlapsSliceRegion = true;
                else
                    hasRemoteSupplementary = true;
            }
        }

        if(!overlapsSliceRegion)
            return;

        if(record.getDuplicateReadFlag())
        {
            if(hasRemoteSupplementary)
                mSliceWriter.registerDuplicateRead(record);
        }
        else if(record.getSupplementaryAlignmentFlag())
        {
            mSliceWriter.addSupplementary(record);
        }
        else
        {
            mSliceWriter.writeRecord(record);
        }
    }

    private boolean overlapsSliceRegion(final String chromosome, final int posStart, final int posEnd)
    {
        if(!mRegion.Chromosome.equals(chromosome))
            return false;

        return mSliceRegions.stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.start(), x.end()));
    }

    private void handleMetricsRecord(final SAMRecord record)
    {

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
        if(mateUnmapped(record))
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
