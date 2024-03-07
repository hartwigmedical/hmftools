package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.samtools.UmiReadType.DUAL;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.common.ReadGroup;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
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
    private final Map<String, ReadGroup> mReadGroupMap; // keyed by readId
    private final CombinedStats mCombinedStats;
    private final ReadCounts mReadCounts;
    private final FlagStats mFlagStats;

    private final List<TargetRegionStats> mTargetRegions;

    private final PerformanceCounter mPerfCounter;
    private final boolean mLogReadIds;
    @Nullable
    private final ChrBaseRegion mPreviousRegion;

    public BamReader(
            final ChrBaseRegion region, final MetricsConfig config, final SamReader samReader, final BamSlicer bamSlicer,
            final CombinedStats combinedStats, @Nullable final ChrBaseRegion previousRegion)
    {
        mConfig = config;
        mRegion = region;
        mPreviousRegion = previousRegion;
        mCombinedStats = combinedStats;

        mSamReader = samReader;
        mBamSlicer = bamSlicer;

        mReadGroupMap = Maps.newHashMap();

        mTargetRegions = Lists.newArrayList();

        if(mConfig.TargetRegions.containsKey(mRegion.Chromosome))
        {
            mConfig.TargetRegions.get(mRegion.Chromosome).stream()
                    .filter(x -> x.overlaps(mRegion))
                    .forEach(x -> mTargetRegions.add(new TargetRegionStats(new ChrBaseRegion(mRegion.Chromosome, x.start(), x.end()))));
        }

        List<ChrBaseRegion> unmappableRegions = mConfig.UnmappableRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

        mBaseCoverage = new BaseCoverage(mConfig, mRegion.start(), mRegion.end(), unmappableRegions);
        mReadCounts = new ReadCounts();
        mFlagStats = new FlagStats();
        mFragmentLengths = new FragmentLengths();

        mPerfCounter = new PerformanceCounter("Slice");
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        BT_LOGGER.debug("processing region({})", mRegion);

        mPerfCounter.start(mConfig.PerfDebug ? mRegion.toString() : null);
        mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);
        mPerfCounter.stop();

        // BT_LOGGER.debug("completed region({})", mRegion);

        postSliceProcess();
    }

    @VisibleForTesting
    protected void postSliceProcess()
    {
        // process overlapping groups
        for(ReadGroup readGroup : mReadGroupMap.values())
        {
            // determine overlapping bases and factor this into the coverage calcs
            processReadGroup(readGroup);
        }

        CoverageMetrics metrics = mBaseCoverage.createMetrics();
        mCombinedStats.addStats(metrics, mFragmentLengths, mReadCounts, mFlagStats, mTargetRegions, mPerfCounter);
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!shouldHandleRead(read))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName()))
        {
            BT_LOGGER.debug("specific readId({}) unmapped({})", read.getReadName(), read.getReadUnmappedFlag());
        }

        processReadForFlagstats(read);
        processReadForFragmentLengths(read);
        processReadForReadCounts(read);
        processReadForTargetRegions(read);
        processReadForBaseCoverage(read);
    }

    private boolean shouldHandleRead(final SAMRecord read)
    {
        if(mRegion.containsPosition(read.getAlignmentStart()))
            return true;

        boolean inPreviousRegion = mPreviousRegion != null && mPreviousRegion.containsPosition(read.getAlignmentStart());
        return !inPreviousRegion && mRegion.overlaps(read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd());
    }

    private void processReadForFlagstats(final SAMRecord read)
    {
        if(isConsensus(read))
        {
            mFlagStats.registerConsensusRead(read);
        }
        else
        {
            mFlagStats.processRead(read);
        }
    }

    private void processReadForFragmentLengths(final SAMRecord read)
    {
        if(!isConsensus(read) && !read.isSecondaryAlignment())
        {
            mFragmentLengths.processRead(read);
        }
    }

    private void processReadForReadCounts(final SAMRecord read)
    {
        addReadToReadCounts(read, mReadCounts);
    }

    private void processReadForTargetRegions(final SAMRecord read)
    {
        for(TargetRegionStats targetRegion : mTargetRegions)
        {
            if(positionsOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), targetRegion.Region.start(), targetRegion.Region.end()))
            {
                addReadToReadCounts(read, targetRegion.Counts);
            }
        }
    }

    private void processReadForBaseCoverage(final SAMRecord read)
    {
        if(read.isSecondaryAlignment())
            return;

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

        if(readHasOverlaps(read))
        {
            readGroup = new ReadGroup(read);
            mReadGroupMap.put(readGroup.id(), readGroup);
            return;
        }

        // process this non-overlapping read immediately without caching
        mBaseCoverage.processRead(read, null);
    }

    private static void addReadToReadCounts(final SAMRecord read, final ReadCounts readCounts)
    {
        if(isConsensus(read))
        {
            // lower the duplicate count to reflect the use of consensus reads
            --readCounts.Duplicates;
            if(isDualStrand(read))
            {
                ++readCounts.DualStrand;
            }
        }
        else
        {
            ++readCounts.TotalReads;
            if(read.getDuplicateReadFlag())
            {
                ++readCounts.Duplicates;
            }
        }
    }

    private static boolean isConsensus(final SAMRecord read)
    {
        return read.hasAttribute(CONSENSUS_READ_ATTRIBUTE);
    }

    private static boolean isDualStrand(final SAMRecord read)
    {
        if(read.getDuplicateReadFlag())
            return false;

        return extractUmiType(read) == DUAL;
    }

    private boolean readHasOverlaps(final SAMRecord read)
    {
        // check the mate and supplementaries for any overlapping bases
        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag() && read.getReferenceName().equals(read.getMateReferenceName()))
        {
            int readLength = read.getReadBases().length;
            int readStart = read.getAlignmentStart();
            int readEnd = read.getAlignmentEnd();

            int mateStart = read.getMateAlignmentStart();
            int mateEnd = mateStart + readLength - 1;

            if(positionsOverlap(readStart, readEnd, mateStart, mateEnd))
                return true;
        }

        if(!read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return false;

        SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

        if(suppData == null)
            return false;

        if(!suppData.Chromosome.equals(read.getReferenceName()))
            return false;

        int suppReadEnd = suppData.Position + read.getReadBases().length - 1;
        return positionsOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), suppData.Position, suppReadEnd);
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
