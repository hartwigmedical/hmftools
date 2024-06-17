package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
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
import com.hartwig.hmftools.common.sv.SvUtils;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

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

    private final List<TargetRegionStats> mTargetRegionStats;
    private final OffTargetFragments mOffTargetFragments;

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

        mTargetRegionStats = Lists.newArrayList();

        if(mConfig.TargetRegions.containsKey(mRegion.Chromosome))
        {
            mConfig.TargetRegions.get(mRegion.Chromosome).stream()
                    .filter(x -> x.overlaps(mRegion))
                    .forEach(x -> mTargetRegionStats.add(new TargetRegionStats(new ChrBaseRegion(mRegion.Chromosome, x.start(), x.end()))));
        }

        mOffTargetFragments = !mConfig.TargetRegions.isEmpty() ? new OffTargetFragments(mConfig.HighFragmentOverlapThreshold) : null;

        List<ChrBaseRegion> unmappableRegions = mConfig.UnmappableRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

        mBaseCoverage = new BaseCoverage(mConfig, mRegion.start(), mRegion.end(), unmappableRegions);
        mReadCounts = new ReadCounts();
        mFlagStats = new FlagStats();
        mFragmentLengths = new FragmentLengths();

        mPerfCounter = new PerformanceCounter("Slice");
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public List<TargetRegionStats> targetRegionStats() { return mTargetRegionStats; }
    public OffTargetFragments offTargetFragments() { return mOffTargetFragments; }

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

            if(read.getDuplicateReadFlag())
                ++mReadCounts.Duplicates;
        }

        mFlagStats.processRead(read, isConsensusRead);

        mFragmentLengths.processRead(read);

        checkTargetRegions(read, isConsensusRead, isDualStrand, readMateEnd);

        if(read.isSecondaryAlignment())
            return;

        // finally handle base coverage (aka depth) - for this any overlaps between reads (including supplementaries) are only counted
        // once and so reads are temporarily cached until their mate is available, and then they can be analysed together
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
            readGroup = new ReadGroup(read, isConsensusRead);
            mReadGroupMap.put(readGroup.id(), readGroup);
            return;
        }

        // process this non-overlapping read immediately without caching
        mBaseCoverage.processRead(read, null, isConsensusRead);
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
            mBaseCoverage.processRead(readGroup.reads().get(0), null, readGroup.IsConsensus);
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
            readGroup.reads().forEach(x -> mBaseCoverage.processRead(x, null, readGroup.IsConsensus));
            return;
        }

        List<int[]> combinedAlignedBaseCoords = Lists.newArrayList();

        for(int i = 0; i < readGroup.reads().size(); ++i)
        {
            SAMRecord read = readGroup.reads().get(i);

            if(i == 0)
            {
                mBaseCoverage.processRead(read, null, readGroup.IsConsensus);
            }
            else
            {
                mBaseCoverage.processRead(read, combinedAlignedBaseCoords, readGroup.IsConsensus);
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
