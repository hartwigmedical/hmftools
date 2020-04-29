package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.ReadRecord;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FusionReadDepth
{
    private final IsofoxConfig mConfig;
    private final Set<String> mChimericReadIds;
    private final Set<String> mReadDepthIds;
    private final FragmentTracker mReadDepthTracker;
    private final List<FusionReadData> mReadDepthFusions;
    private int mReadDepthLimit;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    public FusionReadDepth(final IsofoxConfig config, final Set<String> readIds)
    {
        mConfig = config;
        mChimericReadIds = readIds;
        mReadDepthIds = Sets.newHashSet();
        mReadDepthTracker = new FragmentTracker();
        mReadDepthFusions = Lists.newArrayList();

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = mSamReader != null ? new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true, true) : null;

        mReadDepthLimit = 0;
    }

    private static final int MIN_DEPTH_FRAGS = 1000;
    private static final int MAX_DEPTH_FRAGS = 100000;

    public void calcFusionReadDepth(final List<FusionReadData> fusions)
    {
        // measure the depth at each fusion junction, and at the same time pick up any fragments showing realignment support which
        // where missed the first tme by not being chimeric
        mReadDepthTracker.clear();
        mReadDepthIds.clear();
        mReadDepthFusions.clear();
        mReadDepthFusions.addAll(fusions);

        int maxFusionFrags = fusions.stream().mapToInt(x -> x.fragmentCount()).max().orElse(0);
        mReadDepthLimit = min(MAX_DEPTH_FRAGS, maxFusionFrags * MIN_DEPTH_FRAGS);

        if(mBamSlicer == null)
            return;

        QueryInterval[] queryInterval = new QueryInterval[SE_PAIR];

        int rangeBuffer = mConfig.MaxFragmentLength * 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int[] fusionJuncRange = {-1, -1};

            for(final FusionReadData fusionData : fusions)
            {
                fusionJuncRange[SE_START] = fusionJuncRange[se] == -1 ?
                        fusionData.junctionPositions()[se] : min(fusionJuncRange[SE_START], fusionData.junctionPositions()[se]);

                fusionJuncRange[SE_END] = max(fusionJuncRange[SE_END], fusionData.junctionPositions()[se]);

                // set depth from existing fragments
            }

            int chrSeqIndex = mSamReader.getFileHeader().getSequenceIndex(fusions.get(0).chromosomes()[se]);
            fusionJuncRange[SE_START] -= rangeBuffer;
            fusionJuncRange[SE_END] += rangeBuffer;
            queryInterval[se] = new QueryInterval(chrSeqIndex, fusionJuncRange[SE_START], fusionJuncRange[SE_END]);
        }

        if(queryInterval[0].referenceIndex == queryInterval[1].referenceIndex)
        {
            if(positionsOverlap(queryInterval[0].start, queryInterval[0].end, queryInterval[1].start, queryInterval[1].end))
            {
                QueryInterval[] newQueryInterval = new QueryInterval[1];
                newQueryInterval[0] = new QueryInterval(queryInterval[0].referenceIndex, queryInterval[0].start, queryInterval[1].end);
                queryInterval = newQueryInterval;
            }
            else if(queryInterval[0].start > queryInterval[1].end)
            {
                ISF_LOGGER.error("invalid query regions: lower({}:{} -> {}|) upper({}:{} -> {}|): fusion(count={} first={})",
                        queryInterval[0].referenceIndex, queryInterval[0].start, queryInterval[0].end,
                        queryInterval[0].referenceIndex, queryInterval[1].start, queryInterval[1].end,
                        fusions.size(), fusions.get(0).toString());
                return;
            }
        }

        boolean consumeRecords = false;

        if(consumeRecords)
        {
            mBamSlicer.slice(mSamReader, queryInterval, this::processReadDepthRecord);
        }
        else
        {
            List<SAMRecord> records = mBamSlicer.slice(mSamReader, queryInterval);

            if(records.size() >= 10000)
            {
                final FusionReadData fusionData = mReadDepthFusions.get(0);

                ISF_LOGGER.info("high read depth record count({}) vs maxReads({}), fusion({}) queryRange({}->{} and {}->{})",
                        records.size(), mReadDepthLimit, fusionData.toString(),
                        queryInterval[0].start, queryInterval[0].end,
                        queryInterval.length == 2 ? queryInterval[1].start : queryInterval[0].start,
                        queryInterval.length == 2 ? queryInterval[1].end : queryInterval[0].end);
            }

            for(int i = 0; i < min(mReadDepthLimit, records.size()); ++i)
            {
                processReadDepthRecord(records.get(i));
            }
        }

        // process unpaired fragments for their depth over fusion junctions
        for(Object object : mReadDepthTracker.getValues())
        {
            final ReadRecord read = (ReadRecord) object;

            if(mReadDepthIds.contains(read.Id))
                continue;

            mReadDepthIds.add(read.Id);

            for(final FusionReadData fusionData : mReadDepthFusions)
            {
                checkFragmentDepthSupport(fusionData, Lists.newArrayList(read));
            }
        }

        /*
        for (int se = SE_START; se <= SE_END; ++se)
        {
            specificTest(mReadDepthFusions.get(0).chromosomes()[se], mReadDepthFusions.get(0).junctionPositions()[se]);
        }
        */
    }

    private void processReadDepthRecord(@NotNull final SAMRecord record)
    {
        if(mReadDepthIds.contains(record.getReadName()))
            return; // prevents handling supplementary fragments again

        final ReadRecord read1 = ReadRecord.from(record);
        final ReadRecord read2 = mReadDepthTracker.checkRead(read1);

        if(read2 == null)
            return;

        processReadDepthRecord(read1, read2);

        if(mReadDepthIds.size() > mReadDepthLimit)
        {
            final FusionReadData fusion = mReadDepthFusions.get(0);
            ISF_LOGGER.warn("exiting read depth at limit({}), fusion({}) frags({}) depth({}-{})",
                    mReadDepthIds.size(), fusion.toString(), fusion.getAllFragments().size(),
                    fusion.getReadDepth()[SE_START], fusion.getReadDepth()[SE_END]);
            mBamSlicer.haltProcessing();
            return;
        }
    }

    private void checkFragmentDepthSupport(final FusionReadData fusionData, final List<ReadRecord> reads)
    {
        for (int se = SE_START; se <= SE_END; ++se)
        {
            for (ReadRecord read : reads)
            {
                if (!read.Chromosome.equals(fusionData.chromosomes()[se]))
                    continue;

                final int seIndex = se;
                if (read.getMappedRegionCoords().stream()
                        .anyMatch(x -> positionWithin(fusionData.junctionPositions()[seIndex], x[SE_START], x[SE_END])))
                {
                    ++fusionData.getReadDepth()[se];
                    break;
                }
            }
        }
    }

    @VisibleForTesting
    public void processReadDepthRecord(final ReadRecord read1, final ReadRecord read2)
    {
        mReadDepthIds.add(read1.Id);

        for(final FusionReadData fusionData : mReadDepthFusions)
        {
            checkFragmentDepthSupport(fusionData, Lists.newArrayList(read1, read2));
        }

        if(mChimericReadIds.contains(read1.Id) || !read1.Chromosome.equals(read2.Chromosome))
            return;

        // check for potential realgined fragments which would have been missed initially for not being chimeric
        FusionFragment fragment = new FusionFragment(Lists.newArrayList(read1, read2));

        for(final FusionReadData fusionData : mReadDepthFusions)
        {
            if(fusionData.isRelignedFragment(fragment))
            {
                fragment.setType(REALIGNED);
                fusionData.addFusionFragment(fragment);
            }
        }
    }

    private void specificTest(final String chromosome, int position)
    {
        QueryInterval[] queryInterval = new QueryInterval[1];
        queryInterval[0] = new QueryInterval(mSamReader.getFileHeader().getSequenceIndex(chromosome), position, position);
        List<SAMRecord> records = mBamSlicer.slice(mSamReader, queryInterval);

        Set<String> ids = Sets.newHashSet();

        for(SAMRecord record : records)
        {
            ids.add(record.getReadName());
        }

        for(String readId : ids)
        {
            if(!mReadDepthIds.contains(readId))
            {
                SAMRecord record = records.stream().filter(x -> x.getReadName().equals(readId)).findFirst().orElse(null);
                ReadRecord read = ReadRecord.from(record);
                ISF_LOGGER.debug("missed read depth id({}): ", readId, read);
            }
        }

        /*
        for(String readId : mReadDepthIds)
        {
            if(!ids.contains(readId))
            {
                ISF_LOGGER.debug("extra read depth id({})", readId);
            }
        }

        */
    }

}
