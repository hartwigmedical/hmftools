package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilename;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.writeCombinedFusions;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionData.FILTER_PASS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ISOFOX_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class FusionCohortTask implements Callable
{
    private final int mTaskId;
    private final CohortConfig mConfig;

    private final Map<String,Path> mSampleFileMap;
    private final Map<String,Integer> mFieldsMap;
    private String mFilteredFusionHeader;
    private final FusionFilters mFilters;
    private final BufferedWriter mCombinedFusionWriter;
    private final FusionCollection mFusionCollection;
    private final ExternalFusionCompare mExternalFusionCompare;

    public FusionCohortTask(
            int taskId, final CohortConfig config, final Map<String,Path> sampleFileMap, final FusionFilters filters,
            final FusionCollection fusionCollection, final BufferedWriter combinedFusionWriter, final BufferedWriter extCompareWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mSampleFileMap = sampleFileMap;
        mFieldsMap = Maps.newHashMap();
        mFilteredFusionHeader = null;
        mFilters = filters;
        mFusionCollection = fusionCollection;
        mCombinedFusionWriter = combinedFusionWriter;

        mExternalFusionCompare = extCompareWriter != null ? new ExternalFusionCompare(mConfig, extCompareWriter) : null;
    }

    public final ExternalFusionCompare getExternalCompare() { return mExternalFusionCompare; }

    @Override
    public Long call()
    {
        ISF_LOGGER.info("task {}: processing {} sample fusion files", mTaskId, mSampleFileMap.size());

        int totalProcessed = 0;

        final Map<String,Integer> fieldsMap = Maps.newHashMap();

        // load each sample's fusions and consolidate into a single list
        for(Map.Entry<String,Path> entry : mSampleFileMap.entrySet())
        {
            final String sampleId = entry.getKey();
            final Path fusionFile = entry.getValue();

            ISF_LOGGER.debug("task {}: sample({}:{}) loading fusion data", mTaskId, totalProcessed, sampleId);

            final List<FusionData> sampleFusions = loadSampleFile(fusionFile);

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

            if(!mConfig.Fusions.ComparisonSources.isEmpty())
            {
                mExternalFusionCompare.compareFusions(sampleId, sampleFusions);
            }
        }

        ISF_LOGGER.info("task {}: task complete", mTaskId, mSampleFileMap.size());

        return (long)0;
    }

    private List<FusionData> loadSampleFile(final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
            {
                mFilteredFusionHeader = lines.get(0);
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));
            }

            lines.remove(0);

            List<FusionData> fusions = Lists.newArrayList();

            for(String data : lines)
            {
                FusionData fusion = FusionData.fromCsv(data, mFieldsMap);

                if(mConfig.Fusions.WriteFilteredFusions || mConfig.Fusions.WriteCombinedFusions)
                    fusion.cacheCsvData(data);

                fusions.add(fusion);
            }

            return fusions;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion file({}): {}", filename.toString(), e.toString());
            return Lists.newArrayList();
        }
    }

    private void writeFilteredFusion(final String sampleId, final List<FusionData> sampleFusions)
    {
        // mark passing fusions, and then include any which are related to them
        final List<FusionData> passingFusions = Lists.newArrayList();
        final List<FusionData> nonPassingFusionsWithRelated = Lists.newArrayList();

        for (FusionData fusion : sampleFusions)
        {
            mFilters.markKnownGeneTypes(fusion);

            FusionCohortData cohortMatch = mFilters.findCohortFusion(fusion);

            if(cohortMatch != null)
                fusion.setCohortFrequency(cohortMatch.sampleCount());

            if(mFilters.isPassingFusion(fusion))
            {
                passingFusions.add(fusion);
            }
            else
            {
                if(!fusion.relatedFusionIds().isEmpty())
                    nonPassingFusionsWithRelated.add(fusion);
            }
        }

        int relatedToPassing = 0;
        for (FusionData fusion : nonPassingFusionsWithRelated)
        {
            boolean matchesPassing = false;
            for(FusionData passingFusion : passingFusions)
            {
                if(fusion.isRelated(passingFusion))
                {
                    matchesPassing = true;

                    if(passingFusion.hasKnownSpliceSites())
                    {
                        fusion.setHasRelatedKnownSpliceSites();
                        break;
                    }
                }
            }

            if(matchesPassing)
            {
                ++relatedToPassing;
                passingFusions.add(fusion);
            }
        }

        passingFusions.forEach(x -> x.setFilter(FILTER_PASS));

        ISF_LOGGER.debug("sample({}) passing fusions({}) from total({} relatedToPass={})",
                sampleId, passingFusions.size(), sampleFusions.size(), relatedToPassing);

        writeFusions(sampleId, passingFusions, PASS_FUSION_FILE_ID);

        if(mConfig.Fusions.RewriteAnnotatedFusions)
        {
            writeFusions(sampleId, sampleFusions, FUSION_FILE_ID);
        }

        if(mConfig.Fusions.WriteCombinedFusions)
        {
            writeCombinedFusions(mCombinedFusionWriter, sampleId, passingFusions);
        }
    }

    private void writeFusions(final String sampleId, final List<FusionData> fusions, final String fileId)
    {
        String outputFile = mConfig.OutputDir + sampleId + ISOFOX_ID + fileId;

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write(mFilteredFusionHeader);
            writer.write(",Filter,CohortCount,KnownFusionType");
            writer.newLine();

            for (FusionData fusion : fusions)
            {
                writer.write(fusion.rawData());
                writer.write(String.format(",%s,%d,%s",
                        fusion.filter(), fusion.cohortFrequency(), fusion.getKnownFusionType()));
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
