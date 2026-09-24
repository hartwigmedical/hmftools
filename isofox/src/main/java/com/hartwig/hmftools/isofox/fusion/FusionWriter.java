package com.hartwig.hmftools.isofox.fusion;

import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.WriteType.FUSION_FRAGMENT;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionWriter
{
    private final IsofoxConfig mConfig;
    private BufferedWriter mFusionWriter;
    private BufferedWriter mPassingFusionWriter;
    private BufferedWriter mFragmentWriter;
    private final boolean mWriteFragments;

    private int mNextFusionId;

    public static final String UNFILTERED_FUSION_FILE_ID = "fusions.tsv";

    public FusionWriter(final IsofoxConfig config)
    {
        mConfig = config;
        mWriteFragments = mConfig.WriteTypes.contains(FUSION_FRAGMENT);

        mFusionWriter = null;
        mFragmentWriter = null;
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
    }

    private void initialiseFusionWriters()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            mFusionWriter = createBufferedWriter(mConfig.formOutputFile(UNFILTERED_FUSION_FILE_ID), false);
            mFusionWriter.write(FusionData.header());
            mFusionWriter.newLine();

            mPassingFusionWriter = createBufferedWriter(mConfig.formOutputFile(PASS_FUSION_FILE_ID), false);
            mPassingFusionWriter.write(RnaFusionFile.header());
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
                mFusionWriter.write(fusionData.toTsv());
                mFusionWriter.newLine();
            }

            for(FusionData fusionData : passingFusions)
            {
                RnaFusion fusion = fusionData.buildRnaFusion();
                mPassingFusionWriter.write(RnaFusionFile.write(fusion));
                mPassingFusionWriter.newLine();
            }

            if(mWriteFragments)
            {
                for(List<FusionReadData> fusionCandidate : fusionCandidates.values())
                {
                    for(FusionReadData fusion : fusionCandidate)
                    {
                        for(List<FusionFragment> fragments : fusion.getFragments().values())
                        {
                            for(FusionFragment fragment : fragments)
                            {
                                writeFragmentData(fragment, fusionId(fusion.id()));
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

    public synchronized void writeUnfusedFragments(final List<FusionFragment> fragments)
    {
        if(!mWriteFragments)
            return;

        fragments.forEach(x -> writeFragmentData(x, "UNFUSED"));
    }

    private void initialiseFragmentWriter()
    {
        if(!mWriteFragments)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("fusion_fragment.tsv");

            mFragmentWriter = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("ReadId").add("ReadCount").add("FusionId").add("FragType").add("SameGeneSet").add("ScCount").add("HasSupp");

            // mFragmentWriter.write("ReadId,ReadCount,FusionGroup,Type,SameGeneSet,ScCount,HasSupp");

            for(int se = SE_START; se <= SE_END; ++se)
            {
                String prefix = se == SE_START ? "Start" : "End";
                sj.add("Chr" + prefix);
                sj.add("Orient" + prefix);
                sj.add("JuncPos" + prefix);
                sj.add("JuncOrient" + prefix);
                sj.add("GeneSet" + prefix);
                sj.add("Region" + prefix);
            }

            mFragmentWriter.write(sj.toString());
            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }
    }

    public synchronized void writeFragmentData(final FusionFragment fragment, final String fusionId)
    {
        if(!mWriteFragments)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(fragment.readId());
            sj.add(valueOf(fragment.reads().size()));
            sj.add(fusionId);
            sj.add(valueOf(fragment.type()));
            sj.add(valueOf(fragment.isSingleGeneCollection()));
            sj.add(valueOf(fragment.reads().stream().filter(x -> x.SoftClipLengths[SE_START] > 0 || x.SoftClipLengths[SE_END] > 0).count()));
            sj.add(valueOf(fragment.hasSuppAlignment()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                sj.add(fragment.chromosomes()[se]);
                sj.add(valueOf(fragment.orientations()[se]));
                sj.add(valueOf(fragment.junctionPositions()[se]));
                sj.add(valueOf(fragment.junctionOrientations()[se]));
                sj.add(valueOf(fragment.geneCollections()[se]));
                sj.add(valueOf(fragment.regionMatchTypes()[se]));
            }

            mFragmentWriter.write(sj.toString());
            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write fusion fragment data: {}", e.toString());
        }
    }
}
