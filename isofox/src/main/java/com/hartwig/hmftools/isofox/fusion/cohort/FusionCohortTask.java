package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort.writeCombinedFusions;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.FRAGMENT_COUNT;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.PASS;
import static com.hartwig.hmftools.isofox.fusion.cohort.PassingFusions.hasSufficientKnownFusionFragments;
import static com.hartwig.hmftools.isofox.fusion.cohort.PassingFusions.isShortLocalFusion;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.fusion.FusionData;

public class FusionCohortTask implements Callable
{
    private final int mTaskId;
    private final CohortConfig mConfig;

    private final Map<String,Path> mSampleFileMap;
    private final Map<String,Integer> mFieldsMap;
    private String mFilteredFusionHeader;
    private final PassingFusions mFilters;
    private final BufferedWriter mCombinedFusionWriter;
    private final FusionCollection mFusionCollection;
    private final ExternalFusionCompare mExternalFusionCompare;

    public FusionCohortTask(
            int taskId, final CohortConfig config, final Map<String,Path> sampleFileMap, final PassingFusions filters,
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
    public final Set<String> getSamples() { return mSampleFileMap.keySet(); }

    @Override
    public Long call()
    {
        ISF_LOGGER.info("task {}: processing {} sample fusion files", mTaskId, mSampleFileMap.size());

        int totalProcessed = 0;

        // load each sample's fusions and consolidate into a single list
        for(Map.Entry<String,Path> entry : mSampleFileMap.entrySet())
        {
            final String sampleId = entry.getKey();
            final Path fusionFile = entry.getValue();

            ISF_LOGGER.debug("task {}: sample({}:{}) loading fusion data", mTaskId, totalProcessed, sampleId);

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
        }

        ISF_LOGGER.info("task {}: task complete", mTaskId, mSampleFileMap.size());

        return (long)0;
    }

    /*
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
    */

    private void writeFilteredFusion(final String sampleId, final List<FusionData> sampleFusions)
    {
        // mark passing fusions, and then include any which are related to them
        final List<FusionData> passingFusions = Lists.newArrayList();
        final List<FusionData> nonPassingFusionsWithRelated = Lists.newArrayList();
        final Map<String,List<FusionData>> lowSupportKnownFusions = Maps.newHashMap();

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
                boolean isShortLocal = isShortLocalFusion(fusion);

                if(!isShortLocal && !fusion.relatedFusionIds().isEmpty())
                    nonPassingFusionsWithRelated.add(fusion);

                // combine fragment support for non-local known-pair fusions
                if(fusion.getKnownFusionType() == KnownGeneType.KNOWN_PAIR && fusion.getFilter() == FRAGMENT_COUNT && !isShortLocal)
                {
                    List<FusionData> fusions = lowSupportKnownFusions.get(fusion.name());
                    if(fusions == null)
                        lowSupportKnownFusions.put(fusion.name(), Lists.newArrayList(fusion));
                    else
                        fusions.add(fusion);
                }
            }
        }

        for(final List<FusionData> fusions : lowSupportKnownFusions.values())
        {
            if(hasSufficientKnownFusionFragments(fusions))
            {
                passingFusions.addAll(fusions);
            }
        }

        int relatedToPassing = 0;
        for (FusionData fusion : nonPassingFusionsWithRelated)
        {
            boolean matchesPassing = false;
            for(FusionData passingFusion : passingFusions)
            {
                if(fusion.isRelated(passingFusion) && passingFusion.hasKnownSpliceSites())
                {
                    // use the same filters as for the passing fusion by using it's known and splice types
                    if(mFilters.isPassingFusion(fusion, passingFusion.getKnownFusionType(), passingFusion.hasKnownSpliceSites()))
                    {
                        matchesPassing = true;
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

        passingFusions.forEach(x -> x.setFilter(PASS));

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
        String outputFile = mConfig.OutputDir + sampleId + ISF_FILE_ID + fileId;

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
                        fusion.getFilter(), fusion.cohortFrequency(), fusion.getKnownFusionType()));
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
