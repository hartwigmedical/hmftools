package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionWriter
{
    private final IsofoxConfig mConfig;
    private BufferedWriter mFusionWriter;
    private BufferedWriter mPassingFusionWriter;
    private BufferedWriter mFragmentWriter;
    private final ChimericReadCache mChimericReadCache;
    private final boolean mWriteReads;
    private final boolean mWriteFragments;

    private int mNextFusionId;

    public static final String RAW_FUSION_FILE_ID = "fusions.csv";
    public static final String PASS_FUSION_FILE_ID = "pass_fusions.csv";

    public FusionWriter(final IsofoxConfig config)
    {
        mConfig = config;
        mWriteReads = mConfig.Fusions.WriteChimericReads;
        mWriteFragments = mConfig.Fusions.WriteChimericFragments;

        mFusionWriter = null;
        mFragmentWriter = null;
        mChimericReadCache = new ChimericReadCache(config);
        mNextFusionId = 0;

        initialiseFusionWriters();
        initialiseFragmentWriter();
    }

    public synchronized int getNextFusionId() { return mNextFusionId++; }

    public void close()
    {
        closeBufferedWriter(mFusionWriter);
        closeBufferedWriter(mPassingFusionWriter);
        closeBufferedWriter(mFragmentWriter);
        mChimericReadCache.close();
    }

    private void initialiseFusionWriters()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            mFusionWriter = createBufferedWriter(mConfig.formOutputFile(RAW_FUSION_FILE_ID), false);
            mFusionWriter.write(FusionData.csvHeader(false));
            mFusionWriter.newLine();

            mPassingFusionWriter = createBufferedWriter(mConfig.formOutputFile(PASS_FUSION_FILE_ID), false);
            mPassingFusionWriter.write(FusionData.csvHeader(true));
            mPassingFusionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to create fusions file: {}", e.toString());
        }
    }

    public synchronized void writeFusionData(
            final List<FusionData> fusions, final List<FusionData> passingFusions, final Map<String,List<FusionReadData>> fusionCandidates)
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            for(FusionData fusionData : fusions)
            {
                mFusionWriter.write(fusionData.toCsv(false));
                mFusionWriter.newLine();
            }

            for(FusionData fusionData : passingFusions)
            {
                mPassingFusionWriter.write(fusionData.toCsv(true));
                mPassingFusionWriter.newLine();
            }

            if(mWriteReads || mWriteFragments)
            {
                for(List<FusionReadData> fusionCandidate : fusionCandidates.values())
                {
                    for(FusionReadData fusion : fusionCandidate)
                    {
                        for(List<FusionFragment> fragments : fusion.getFragments().values())
                        {
                            for(FusionFragment fragment : fragments)
                            {
                                if(mWriteFragments)
                                    writeFragmentData(fragment, fusionId(fusion.id()));

                                if(mWriteReads)
                                    writeReadData(fragment.readId(), fragment.reads(), fusionId(fusion.id()));
                            }
                        }
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    public synchronized void writeReadData(final String readId, final List<FusionRead> reads, final String groupStatus)
    {
        if(mWriteReads)
            mChimericReadCache.writeReadData(readId, reads, groupStatus);
    }

    public synchronized void writeUnfusedFragments(final List<FusionFragment> fragments)
    {
        if(!mWriteFragments)
            return;

        fragments.forEach(x -> writeFragmentData(x, "UNFUSED"));
        fragments.forEach(x -> writeReadData(x.readId(), x.reads(), "UNFUSED"));
    }

    private void initialiseFragmentWriter()
    {
        if(!mWriteFragments)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("chimeric_frags.csv");

            mFragmentWriter = createBufferedWriter(outputFileName, false);
            mFragmentWriter.write("ReadId,ReadCount,FusionGroup,Type,SameGeneSet,ScCount,HasSupp");

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final String prefix = se == SE_START ? "Start" : "End";
                mFragmentWriter.write(",Chr" + prefix);
                mFragmentWriter.write(",Orient" + prefix);
                mFragmentWriter.write(",JuncPos" + prefix);
                mFragmentWriter.write(",JuncOrient" + prefix);
                mFragmentWriter.write(",JuncType" + prefix);
                mFragmentWriter.write(",GeneSet" + prefix);
                mFragmentWriter.write(",Region" + prefix);
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }
    }

    public void writeIncompleteGroupReads(final List<FusionReadGroup> incompleteGroups)
    {
        incompleteGroups.forEach(x -> mChimericReadCache.writeReadData(x.ReadId, x.Reads, "INCOMPLETE_GROUPS"));
    }

    public synchronized void writeFragmentData(final FusionFragment fragment, final String fusionId)
    {
        if(!mWriteFragments)
            return;

        try
        {
            mFragmentWriter.write(String.format("%s,%d,%s,%s,%s,%d,%s",
                    fragment.readId(), fragment.reads().size(), fusionId, fragment.type(),
                    fragment.isSingleGeneCollection(),
                    fragment.reads().stream().filter(x -> x.SoftClipLengths[SE_START] > 0 || x.SoftClipLengths[SE_END] > 0).count(),
                    fragment.hasSuppAlignment()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                mFragmentWriter.write(String.format(",%s,%d,%d,%d,%s,%d,%s",
                        fragment.chromosomes()[se], fragment.orientations()[se],
                        fragment.junctionPositions()[se], fragment.junctionOrientations()[se], fragment.junctionTypes()[se],
                        fragment.geneCollections()[se], fragment.regionMatchTypes()[se]));
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }

    }

}
