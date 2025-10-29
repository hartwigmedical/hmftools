package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.rna.GeneFusionFile.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.common.rna.GeneFusionFile.UNFILTERED_FUSION_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.writeCombinedFusions;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.fusion.FusionData;
import com.hartwig.hmftools.isofox.fusion.PassingFusions;

public class FusionCohortTask implements Callable<Void>
{
    private final int mTaskId;
    private final CohortConfig mConfig;

    private final Map<String, Path> mSampleFileMap;
    private final PassingFusions mFilters;
    private final BufferedWriter mCombinedFusionWriter;
    private final FusionCollection mFusionCollection;
    private final ExternalFusionCompare mExternalFusionCompare;
    private final UnknownSpliceAnalyser mUnknownSpliceAnalyser;

    public FusionCohortTask(
            int taskId, final CohortConfig config, final Map<String,Path> sampleFileMap, final PassingFusions filters,
            final FusionCollection fusionCollection, final BufferedWriter combinedFusionWriter, final BufferedWriter extCompareWriter,
            final BufferedWriter unknownSpliceWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mSampleFileMap = sampleFileMap;
        mFilters = filters;
        mFusionCollection = fusionCollection;
        mCombinedFusionWriter = combinedFusionWriter;

        mExternalFusionCompare = extCompareWriter != null ? new ExternalFusionCompare(mConfig, extCompareWriter) : null;
        mUnknownSpliceAnalyser = unknownSpliceWriter != null ? new UnknownSpliceAnalyser(mConfig.Fusions, unknownSpliceWriter) : null;
    }

    public final ExternalFusionCompare getExternalCompare() { return mExternalFusionCompare; }
    public final Set<String> getSamples() { return mSampleFileMap.keySet(); }

    @Override
    public Void call()
    {
        ISF_LOGGER.info("task {}: processing {} sample fusion files", mTaskId, mSampleFileMap.size());

        int totalProcessed = 0;

        // load each sample's fusions and consolidate into a single list
        for(Map.Entry<String, Path> entry : mSampleFileMap.entrySet())
        {
            final String sampleId = entry.getKey();
            final Path fusionFile = entry.getValue();

            // ISF_LOGGER.debug("task {}: sample({}:{}) loading fusion data", mTaskId, totalProcessed, sampleId);

            final List<FusionData> sampleFusions = FusionData.loadFromFile(fusionFile);

            ISF_LOGGER.info("task {}: sample({}:{}) loaded {} fusions", mTaskId, totalProcessed, sampleId, sampleFusions.size());

            ++totalProcessed;

            if(mConfig.Fusions.WriteFilteredFusions)
            {
                writeFilteredFusion(sampleId, sampleFusions);
                continue;
            }

            // add to the fusion cohort collection
            if(mConfig.Fusions.GenerateCohort)
            {
                mFusionCollection.addToCohortCache(sampleFusions, sampleId);
            }

            if(mConfig.Fusions.ComparisonSource != null)
            {
                mExternalFusionCompare.compareFusions(sampleId, sampleFusions);
            }

            if(mConfig.Fusions.FindUnknownSplice)
            {
                mUnknownSpliceAnalyser.compareFusions(sampleId, sampleFusions);
            }
        }

        ISF_LOGGER.info("task {}: task complete", mTaskId, mSampleFileMap.size());

        return null;
    }

    private void writeFilteredFusion(final String sampleId, final List<FusionData> sampleFusions)
    {
        // mark passing fusions, and then include any which are related to them
        final List<FusionData> passingFusions = mFilters.findPassingFusions(sampleFusions);

        ISF_LOGGER.debug("sample({}) passing fusions({}) from total({})",
                sampleId, passingFusions.size(), sampleFusions.size());

        writeFusions(sampleId, passingFusions, true);

        if(mConfig.Fusions.RewriteAnnotatedFusions)
        {
            writeFusions(sampleId, sampleFusions, false);
        }

        if(mConfig.Fusions.WriteCombinedFusions)
        {
            writeCombinedFusions(mCombinedFusionWriter, sampleId, passingFusions);
        }
    }

    private void writeFusions(final String sampleId, final List<FusionData> fusions, boolean includeFilterFields)
    {
        final String fileId = includeFilterFields ? PASS_FUSION_FILE_ID : UNFILTERED_FUSION_FILE_ID;
        String outputFile = mConfig.OutputDir + sampleId + ISF_FILE_ID + fileId;

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write(FusionData.csvHeader(includeFilterFields));
            writer.newLine();

            for(FusionData fusion : fusions)
            {
                writer.write(fusion.toCsv(includeFilterFields));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write filtered fusion file({}): {}", outputFile, e.toString());
        }
    }
}
