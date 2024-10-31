package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.bam.UmiReadType.DUAL;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.sv.SvUtils;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class BamReader
{
    private final MetricsConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final BaseCoverage mBaseCoverage;
    private final FragmentLengths mFragmentLengths;
    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final CombinedStats mCombinedStats;
    private final ReadCounts mReadCounts;
    private final FlagStats mFlagStats;
    private final PartitionStats mPartitionStats;

    private final List<TargetRegionStats> mTargetRegionStats;
    private final OffTargetFragments mOffTargetFragments;
    private final ChrBaseRegion mExcludedRegion;

    private final PerformanceCounter mPerfCounter;
    private boolean mLogReadIds;
    private long mLastFragmentMapReadCountCheck;

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

        mTargetRegionStats = Lists.newArrayList();

        if(mConfig.TargetRegions.containsKey(mRegion.Chromosome))
        {
            mConfig.TargetRegions.get(mRegion.Chromosome).stream()
                    .filter(x -> x.overlaps(mRegion))
                    .forEach(x -> mTargetRegionStats.add(new TargetRegionStats(new ChrBaseRegion(mRegion.Chromosome, x.start(), x.end()))));
        }

        mOffTargetFragments = !mConfig.TargetRegions.isEmpty() ? new OffTargetFragments(mConfig.HighFragmentOverlapThreshold) : null;

        List<ChrBaseRegion> unmappableRegions = mConfig.UnmappableRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

        ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(mConfig.RefGenVersion);
        mExcludedRegion = excludedRegion.overlaps(mRegion) ? excludedRegion : null;

        mBaseCoverage = new BaseCoverage(mConfig, mRegion.start(), mRegion.end(), unmappableRegions);
        mReadCounts = new ReadCounts();
        mFlagStats = new FlagStats();
        mFragmentLengths = new FragmentLengths();
        mPartitionStats = new PartitionStats();

        mPerfCounter = new PerformanceCounter("Slice");
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mLastFragmentMapReadCountCheck = 0;
    }

    public ChrBaseRegion region() { return mRegion; }
    public List<TargetRegionStats> targetRegionStats() { return mTargetRegionStats; }
    public OffTargetFragments offTargetFragments() { return mOffTargetFragments; }
    public PartitionStats partitionStats() { return mPartitionStats; }

    public void run()
    {
        BT_LOGGER.debug("processing region({})", mRegion);

        mPerfCounter.start(mConfig.PerfDebug ? mRegion.toString() : null);
        mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);
        mPerfCounter.stop();

        mPartitionStats.ProcessTime = mPerfCounter.getLastTime();

        // BT_LOGGER.debug("completed region({})", mRegion);

        postSliceProcess();
    }

    @VisibleForTesting
    protected void postSliceProcess()
    {
        mReadGroupMap.clear();

        CoverageMetrics metrics = mBaseCoverage.createMetrics();

        mCombinedStats.addStats(
                metrics, mFragmentLengths, mReadCounts, mFlagStats,
                mOffTargetFragments != null ? mOffTargetFragments.fragmentOverlapCounts() : Collections.EMPTY_MAP, mPerfCounter);
    }

    private void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();

        if(!mRegion.containsPosition(readStart) && !mConfig.OnlyTargetRegions)
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName()))
        {
            BT_LOGGER.debug("specific readId({}) unmapped({})", read.getReadName(), read.getReadUnmappedFlag());
        }

        // extract mate alignment since used by a few components
        int readMateEnd = read.getReadPairedFlag() ? SamRecordUtils.getMateAlignmentEnd(read) : NO_POSITION;

        boolean isDualStrand = false;

        boolean isConsensusRead = read.hasAttribute(CONSENSUS_READ_ATTRIBUTE);

        if(isConsensusRead)
        {
            isDualStrand = extractUmiType(read) == DUAL;

            // lower the duplicate count to reflect the use of consensus reads - since 2 duplicates would usually result in one being
            // marked as primary and the other as duplicate, whereas now both are duplicates plus a new consensus read are expected
            --mReadCounts.Duplicates;

            if(isDualStrand)
                ++mReadCounts.DualStrand;
        }
        else
        {
            ++mReadCounts.Total;
            checkLogPartitionReadCount();
            purgeReadGroups(readStart);

            if(read.getDuplicateReadFlag())
                ++mReadCounts.Duplicates;
        }

        mFlagStats.processRead(read, isConsensusRead);
        mFragmentLengths.processRead(read);
        mPartitionStats.processRead(read, mRegion, isConsensusRead);

        checkTargetRegions(read, isConsensusRead, isDualStrand, readMateEnd);

        if(read.isSecondaryAlignment())
            return;

        boolean inExcludedRegion = mExcludedRegion != null && positionsOverlap(
                mExcludedRegion.start(), mExcludedRegion.end(), readStart, read.getAlignmentEnd());

        if(inExcludedRegion)
            return;

        // finally handle base coverage (aka depth) - for this any overlaps between reads (including supplementaries) are only counted
        // once and so read coordinates are temporarily cached until any further overlapping reads for the same fragment have been processed
        ReadGroup readGroup = mReadGroupMap.get(read.getReadName());

        if(readGroup != null && readStart >= readGroup.MaxReadStart)
        {
            mReadGroupMap.remove(read.getReadName());
        }

        updateBaseCoverage(read, isConsensusRead, readGroup);
    }

    private void checkTargetRegions(final SAMRecord read, boolean isConsensus, boolean isDualStrand, int readMateEnd)
    {
        if(mTargetRegionStats.isEmpty() && mOffTargetFragments == null)
            return;

        if(read.getSupplementaryAlignmentFlag())
            return;

        boolean isLocalConcordantFragment;

        // note that local concordant fragments are counted once, so totals are fragment totals but for discordant fragments
        // in which case the read is counted twice, ie remotely too as on or off target
        if(read.getReadPairedFlag())
        {
            isLocalConcordantFragment = !SvUtils.isDiscordant(read);

            if(isLocalConcordantFragment)
            {
                // avoid double-counting
                if(read.getAlignmentStart() > read.getMateAlignmentStart())
                    return;
                else if(read.getAlignmentStart() == read.getMateAlignmentStart() && read.getSecondOfPairFlag())
                    return;
            }
        }
        else
        {
            isLocalConcordantFragment = false;
        }

        int alignedPosStart = isLocalConcordantFragment ? min(read.getAlignmentStart(), read.getMateAlignmentStart()) : read.getAlignmentStart();
        int alignedPosEnd = isLocalConcordantFragment ? max(read.getAlignmentEnd(), readMateEnd) : read.getAlignmentEnd();

        boolean overlapsTargetRegion = false;

        for(TargetRegionStats targetRegion : mTargetRegionStats)
        {
            if(!positionsOverlap(alignedPosStart, alignedPosEnd, targetRegion.start(), targetRegion.end()))
                continue;

            overlapsTargetRegion = true;

            if(isConsensus)
            {
                --targetRegion.FragmentCounts.Duplicates;
            }
            else
            {
                ++targetRegion.FragmentCounts.Total;
            }

            if(read.getDuplicateReadFlag())
            {
                ++targetRegion.FragmentCounts.Duplicates;
            }
            else
            {
                if(isDualStrand)
                   ++targetRegion.FragmentCounts.DualStrand;
            }
        }

        if(!read.getDuplicateReadFlag() && !overlapsTargetRegion && mOffTargetFragments != null)
            mOffTargetFragments.addRead(read, isLocalConcordantFragment, readMateEnd);
    }

    private int findOverlappingFragmentStart(final SAMRecord read)
    {
        // check the mate and supplementaries for any overlapping bases - cache read if they are still expected
        int maxGroupReadStart = 0;
        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag() && read.getReferenceName().equals(read.getMateReferenceName()))
        {
            int mateStart = read.getMateAlignmentStart();

            if(mateStart >= readStart && positionWithin(mateStart, readStart, readEnd))
                maxGroupReadStart = mateStart;
        }

        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

            if(suppData != null && suppData.Chromosome.equals(read.getReferenceName()))
            {
                if(suppData.Position >= suppData.Position && positionWithin(suppData.Position, readStart, readEnd))
                    maxGroupReadStart = max(suppData.Position, maxGroupReadStart);
            }
        }

        if(maxGroupReadStart > 0 && mRegion.containsPosition(maxGroupReadStart))
            return maxGroupReadStart;

        return -1;
    }

    private void updateBaseCoverage(final SAMRecord read, final boolean isConsensus, @Nullable final ReadGroup readGroup)
    {
        List<int[]> alignedBaseCoords = readGroup != null ? readGroup.CombinedAlignedBaseCoords : null;

        mBaseCoverage.processRead(read, alignedBaseCoords, isConsensus);

        if(readGroup == null)
        {
            int maxGroupReadStart = findOverlappingFragmentStart(read);

            if(maxGroupReadStart <= 0)
                return;

            ReadGroup newReadGroup = new ReadGroup(read, maxGroupReadStart, isConsensus);
            mReadGroupMap.put(read.getReadName(), newReadGroup);
            alignedBaseCoords = newReadGroup.CombinedAlignedBaseCoords;
        }

        addAlignedCoords(read, alignedBaseCoords);
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
                    int coordsStart = position;
                    int coordsEnd = position + element.getLength() - 1;

                    if(alignedBaseCoords.stream().noneMatch(x -> x[0] == coordsStart && x[1] == coordsEnd))
                        alignedBaseCoords.add(new int[] { coordsStart, coordsEnd});

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

    private void checkLogPartitionReadCount()
    {
        if(mConfig.PartitionReadCountCheck == 0)
            return;

        if(mPartitionStats.TotalReads > 0 && (mPartitionStats.TotalReads % mConfig.PartitionReadCountCheck) == 0)
        {
            BT_LOGGER.debug("partition({}) reads(total={} chimeric={}, interpartition={}) cachedFragments({})",
                    mRegion, mPartitionStats.TotalReads, mPartitionStats.ChimericReads, mPartitionStats.InterPartition,
                    mReadGroupMap.size());
        }
    }

    private static final int FRAGMENT_MAP_SIZE_READ_COUNT_CHECK = 10000;
    private static final int FRAGMENT_MAP_SIZE_THRESHOLD = 10000;

    private void purgeReadGroups(final int currentReadStartPosition)
    {
        if(mReadGroupMap.size() < FRAGMENT_MAP_SIZE_THRESHOLD)
            return;

        if(mPartitionStats.TotalReads < mLastFragmentMapReadCountCheck + FRAGMENT_MAP_SIZE_READ_COUNT_CHECK)
            return;

        mLastFragmentMapReadCountCheck = mPartitionStats.TotalReads;

        // a fall-back routine - remove
        List<String> groupsToRemove = mReadGroupMap.entrySet().stream()
                .filter(x -> x.getValue().MaxReadStart < currentReadStartPosition)
                .map(x -> x.getKey()).collect(Collectors.toList());

        if(groupsToRemove.isEmpty())
            return;

        groupsToRemove.forEach(x -> mReadGroupMap.remove(x));

        BT_LOGGER.trace("partition({}) purged {} readGroups, cachedFragments({})",
                mRegion, groupsToRemove.size(), mReadGroupMap.size());
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
