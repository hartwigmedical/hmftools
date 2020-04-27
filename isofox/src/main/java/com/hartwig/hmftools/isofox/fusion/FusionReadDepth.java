package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FusionReadDepth
{
    private final IsofoxConfig mConfig;
    private final Map<String,List<ReadRecord>> mReadsMap;
    private final FragmentTracker mReadDepthTracker;
    private final List<FusionReadData> mReadDepthFusions;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    public FusionReadDepth(final IsofoxConfig config, final Map<String,List<ReadRecord>> readsMap)
    {
        mConfig = config;
        mReadsMap = readsMap;
        mReadDepthTracker = new FragmentTracker();
        mReadDepthFusions = Lists.newArrayList();

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = mSamReader != null ? new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true, true) : null;
    }

    public void calcFusionReadDepth(final List<FusionReadData> fusions)
    {
        mReadDepthTracker.clear();
        mReadDepthFusions.clear();
        mReadDepthFusions.addAll(fusions);

        if(mBamSlicer == null)
            return;

        QueryInterval[] queryInterval = new QueryInterval[SE_PAIR];

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
            queryInterval[se] = new QueryInterval(chrSeqIndex, fusionJuncRange[SE_START], fusionJuncRange[SE_END]);
        }

        if(queryInterval[0].referenceIndex == queryInterval[1].referenceIndex && queryInterval[0].start == queryInterval[1].start
        && queryInterval[0].end == queryInterval[1].end)
        {
            ISF_LOGGER.error("invalid query region: fusion(count={} first={})", fusions.size(), fusions.get(0).toString());
            return;
        }

        mBamSlicer.slice(mSamReader, queryInterval, this::processReadDepthRecord);

        // process unpaired fragments for their depth over fusion junctions
        for(Object object : mReadDepthTracker.getValues())
        {
            final ReadRecord read = (ReadRecord) object;

            for(final FusionReadData fusionData : mReadDepthFusions)
            {
                for (int se = SE_START; se <= SE_END; ++se)
                {
                    final List<int[]> readCoords = read.getMappedRegionCoords();

                    if (!read.Chromosome.equals(fusionData.chromosomes()[se]))
                        continue;

                    final int seIndex = se;

                    if (readCoords.stream().anyMatch(x -> positionWithin(fusionData.junctionPositions()[seIndex], x[SE_START], x[SE_END])))
                    {
                        ++fusionData.getReadDepth()[se];
                    }
                }
            }
        }
    }

    private void processReadDepthRecord(@NotNull final SAMRecord record)
    {
        if(mReadsMap.containsKey(record.getReadName()))
            return;

        final ReadRecord read1 = ReadRecord.from(record);
        final ReadRecord read2 = mReadDepthTracker.checkRead(read1);

        if(read2 == null)
            return;

        processReadDepthRecord(read1, read2);
    }

    @VisibleForTesting
    public void processReadDepthRecord(final ReadRecord read1, final ReadRecord read2)
    {
        if(!read1.Chromosome.equals(read2.Chromosome))
        {
            ISF_LOGGER.warn("read-depth read({}) spans chromosomes({} & {})", read1.Id, read1.Chromosome, read2.Chromosome);
            return;
        }

        // test these pairs against each fusion junction
        FusionFragment fragment = new FusionFragment(Lists.newArrayList(read1, read2));

        for(final FusionReadData fusionData : mReadDepthFusions)
        {
            if(fusionData.isRelignedFragment(fragment))
            {
                fragment.setType(REALIGNED);
                fusionData.addFusionFragment(fragment);
                continue;
            }

            boolean[] hasReadSupport = { false, false };

            for (int se = SE_START; se <= SE_END; ++se)
            {
                for (ReadRecord read : fragment.getReads())
                {
                    final List<int[]> readCoords = read.getMappedRegionCoords();

                    if (!read.Chromosome.equals(fusionData.chromosomes()[se]))
                        continue;

                    final int seIndex = se;

                    if (!hasReadSupport[se] && readCoords.stream()
                            .anyMatch(x -> positionWithin(fusionData.junctionPositions()[seIndex], x[SE_START], x[SE_END])))
                    {
                        hasReadSupport[se] = true;
                    }
                }
            }

            for (int se = SE_START; se <= SE_END; ++se)
            {
                if (hasReadSupport[se])
                    ++fusionData.getReadDepth()[se];
            }
        }
    }
}
