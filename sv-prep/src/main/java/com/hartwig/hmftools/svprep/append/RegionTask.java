package com.hartwig.hmftools.svprep.append;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvCommon.createBamSlicer;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.append.AppendConstants.BREAKEND_PROXIMITY;
import static com.hartwig.hmftools.svprep.append.AppendConstants.JUNCTION_DISTANCE_BUFFER;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;
import com.hartwig.hmftools.svprep.HotspotCache;
import com.hartwig.hmftools.svprep.reads.JunctionData;
import com.hartwig.hmftools.svprep.reads.JunctionTracker;
import com.hartwig.hmftools.svprep.reads.JunctionsConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilterType;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionTask implements Callable
{
    private final AppendConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final List<BreakendData> mBreakends;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    // state
    private int mTotalReads;
    private final JunctionTracker mBreakendTracker;

    private boolean mLogReadIds;

    public RegionTask(final AppendConfig config, final String chromosome, final List<BreakendData> breakends)
    {
        mConfig = config;
        mBreakends = Lists.newArrayList(breakends);

        int regionStart = mBreakends.get(0).Position - BREAKEND_PROXIMITY;
        int regionEnd = mBreakends.get(mBreakends.size() - 1).Position + BREAKEND_PROXIMITY;

        mRegion = new ChrBaseRegion(chromosome, regionStart, regionEnd);

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
        mBamSlicer = createBamSlicer();

        JunctionsConfig jtConfig = new JunctionsConfig(
                true, mConfig.ReadFiltering, DEFAULT_READ_LENGTH, false, false,
                false, true, 0, false);

        mBreakendTracker = new JunctionTracker(mRegion, jtConfig, new HotspotCache(null), new BlacklistLocations(null));

        // register existing breakends as junctions
        List<JunctionData> junctions = breakends.stream()
                .map(x -> new JunctionData(x.Position, x.Orientation, null)).collect(Collectors.toList());

        mBreakendTracker.addExistingJunctions(junctions);

        mTotalReads = 0;

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    @Override
    public Long call()
    {
        SV_LOGGER.debug("slicing region({}) for {} breakends", mRegion, mBreakends.size());

        mBamSlicer.slice(mSamReader, mRegion, this::processRead);

        mBreakendTracker.assignFragments();

        SV_LOGGER.debug("region({}) complete, total reads({})", mRegion, mTotalReads);

        for(BreakendData breakendData : mBreakends)
        {
            List<JunctionData> junctions = mBreakendTracker.junctions().stream()
                    .filter(x -> abs(x.Position - breakendData.Position) <= JUNCTION_DISTANCE_BUFFER
                            && x.Orientation == breakendData.Orientation)
                    .collect(Collectors.toList());

            breakendData.addJunctionData(junctions);
        }

        mBreakendTracker.clear();

        return (long)0;
    }

    private void processRead(final SAMRecord record)
    {
        ++mTotalReads;

        int readStart = record.getAlignmentStart();

        if(readStart > mRegion.end())
        {
            mBamSlicer.haltProcessing();
            return;
        }

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(record.getReadName()))
                SV_LOGGER.debug("specific readId({})", record.getReadName());
        }

        int filters = mConfig.ReadFiltering.checkFilters(record);

        if(filters == 0 || filters == ReadFilterType.MIN_MAP_QUAL.flag())
        {
            ReadRecord read = ReadRecord.from(record);
            read.setFilters(filters);
            read.setReadType(JUNCTION);

            mBreakendTracker.processRead(read);
        }
        else
        {
            boolean isSupportCandidate = mConfig.ReadFiltering.isCandidateSupportingRead(record, filters);

            if(!isSupportCandidate)
                return;

            ReadRecord read = ReadRecord.from(record);
            read.setFilters(filters);

            if(isSupportCandidate)
                read.setReadType(CANDIDATE_SUPPORT);

            mBreakendTracker.processRead(read);
        }
    }

}
