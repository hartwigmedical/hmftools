package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.ONE_JUNCTION;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
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
    private final Map<String,List<ReadRecord>> mReadsMap;
    private final FragmentTracker mReadDepthTracker;
    private final List<FusionReadData> mReadDepthFusions;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final static int SOFT_CLIP_JUNCTION_DISTANCE = 5;

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
        if(mBamSlicer == null)
            return;

        mReadDepthTracker.clear();
        mReadDepthFusions.clear();
        mReadDepthFusions.addAll(fusions);

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

        mBamSlicer.slice(mSamReader, queryInterval, this::processReadDepthRecord);
    }

    private void processReadDepthRecord(@NotNull final SAMRecord record)
    {
        if(mReadsMap.containsKey(record.getReadName()))
            return;

        final ReadRecord read1 = ReadRecord.from(record);
        final ReadRecord read2 = mReadDepthTracker.checkRead(read1);

        if(read2 == null)
            return;

        if(!read1.containsSoftClipping() && !read2.containsSoftClipping())
            return;

        if(!read1.Chromosome.equals(read2.Chromosome))
        {
            ISF_LOGGER.warn("read-depth read({}) spans chromosomes({} & {})", read1.Id, read1.Chromosome, read2.Chromosome);
            return;
        }

        // test these pairs against each fusion junction
        FusionFragment fragment = new FusionFragment(Lists.newArrayList(read1, read2));

        for(final FusionReadData fusionData : mReadDepthFusions)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean hasReadSupport = false;
                boolean addedFragment = false;

                for(ReadRecord read : fragment.getReads())
                {
                    final List<int[]> readCoords = read.getMappedRegionCoords();
                    final int[] readBoundaries = { readCoords.get(0)[SE_START], readCoords.get(readCoords.size() - 1)[SE_END] };

                    if(!read.Chromosome.equals(fusionData.chromosomes()[se]))
                        continue;

                    final int seIndex = se;

                    if(!hasReadSupport && readCoords.stream()
                            .anyMatch(x -> positionWithin(fusionData.junctionPositions()[seIndex], x[SE_START], x[SE_END])))
                    {
                        hasReadSupport = true;
                    }

                    if(fusionData.junctionOrientations()[se] == 1 && read.isSoftClipped(SE_END))
                    {
                        if(positionWithin(
                                fusionData.junctionPositions()[se],
                                readBoundaries[SE_END], readBoundaries[SE_END] + SOFT_CLIP_JUNCTION_DISTANCE))
                        {
                            fragment.setType(ONE_JUNCTION);
                            fragment.setGeneData(fusionData.geneCollections()[se], fusionData.getTransExonRefsByPos(se));
                            fusionData.addFusionFragment(fragment);
                            addedFragment = true;
                            break;
                        }
                    }
                    else if(fusionData.junctionOrientations()[se] == -1 && read.isSoftClipped(SE_START))
                    {
                        if(positionWithin(
                                fusionData.junctionPositions()[se],
                                readBoundaries[SE_START] - SOFT_CLIP_JUNCTION_DISTANCE, readBoundaries[SE_START]))
                        {
                            fragment.setType(ONE_JUNCTION);
                            fragment.setGeneData(fusionData.geneCollections()[se], fusionData.getTransExonRefsByPos(se));
                            fusionData.addFusionFragment(fragment);
                            addedFragment = true;
                            break;
                        }
                    }
                }

                if(hasReadSupport && !addedFragment)
                    ++fusionData.getReadDepth()[se];
            }
        }
    }
}
