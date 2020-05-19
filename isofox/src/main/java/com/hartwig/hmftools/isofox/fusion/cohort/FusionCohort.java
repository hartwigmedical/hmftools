package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.PASSING_FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortAnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;

public class FusionCohort
{
    private final CohortConfig mConfig;

    private final FusionFilters mFilters;
    private final Map<String,Integer> mFieldsMap;
    private String mFilteredFusionHeader;
    private final ExternalFusionCompare mExternalFusionCompare;

    // map of chromosome-pair to start position to list of fusion junctions
    private final Map<String, Map<Integer,List<FusionCohortData>>> mFusions;

    private int mFusionCount;

    public static final String PASS_FUSION_FILE_ID = "pass_fusions.csv";

    /* Routines:
        1. generate a cohort file from multiple sample fusion files from Isofox
        2. write a set of filtered/passing fusions for each sample fusion file loaded
        3. run comparison of filtered/passing fusions between Isofox and external fusion files
    */

    public FusionCohort(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mFilters = new FusionFilters(mConfig.Fusions, cmd);
        mFusions = Maps.newHashMap();
        mFieldsMap = Maps.newHashMap();
        mFusionCount = 0;
        mFilteredFusionHeader = null;
        mExternalFusionCompare = !mConfig.Fusions.ComparisonSources.isEmpty() ? new ExternalFusionCompare(mConfig) : null;
    }

    public void processFusionFiles()
    {
        if(!mConfig.Fusions.GenerateCohort && mConfig.Fusions.ComparisonSources.isEmpty() && !mConfig.Fusions.WriteFilteredFusions)
        {
            ISF_LOGGER.warn("no fusion functions configured");
            return;
        }

        final List<Path> filenames = Lists.newArrayList();

        final CohortAnalysisType fileType = !mConfig.Fusions.ComparisonSources.isEmpty() ? PASSING_FUSION : FUSION;

        if(!formSampleFilenames(mConfig, fileType, filenames))
            return;

        ISF_LOGGER.info("loading {} sample fusion files", mConfig.SampleData.SampleIds.size());

        int totalProcessed = 0;
        int nextLog = 100000;

        // load each sample's fusions and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path fusionFile = filenames.get(i);

            ISF_LOGGER.debug("{}: sample({}) loading fusion data", i, sampleId);

            final List<FusionData> sampleFusions = loadSampleFile(fusionFile);

            ISF_LOGGER.info("{}: sample({}) loaded {} fusions", i, sampleId, sampleFusions.size());

            if(mConfig.Fusions.WriteFilteredFusions)
            {
                writeFilteredFusion(sampleId, sampleFusions);
                continue;
            }

            // add to the fusion cohort collection
            if(mConfig.Fusions.GenerateCohort)
            {
                sampleFusions.forEach(x -> addToCohortCache(x, sampleId));
            }

            if(!mConfig.Fusions.ComparisonSources.isEmpty())
            {
                mExternalFusionCompare.compareFusions(sampleId, sampleFusions);
            }

            if(mFusionCount > nextLog)
            {
                nextLog += 100000;
                ISF_LOGGER.info("total fusion count({})", mFusionCount);
            }
        }

        if(mConfig.Fusions.GenerateCohort)
        {
            ISF_LOGGER.info("loaded {} fusion records, total fusion count({})", totalProcessed, mFusionCount);
            FusionCohortData.writeCohortFusions(mFusions, mConfig, mConfig.Fusions);
            return;
        }

        if(mExternalFusionCompare != null)
            mExternalFusionCompare.close();
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

                if(mConfig.Fusions.WriteFilteredFusions)
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

    private void addToCohortCache(final FusionData fusion, final String sampleId)
    {
        if(fusion.SplitFrags < mConfig.Fusions.MinFragCount)
            return;

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            if(!mConfig.RestrictedGeneIds.contains(fusion.GeneIds[SE_START]) || !mConfig.RestrictedGeneIds.contains(fusion.GeneIds[SE_END]))
                return;
        }

        if(!mConfig.ExcludedGeneIds.isEmpty())
        {
            if(mConfig.ExcludedGeneIds.contains(fusion.GeneIds[SE_START]) || mConfig.ExcludedGeneIds.contains(fusion.GeneIds[SE_END]))
                return;
        }

        final String chrPair = formChromosomePair(fusion.Chromosomes[SE_START], fusion.Chromosomes[SE_END]);

        Map<Integer,List<FusionCohortData>> chrPairFusions = mFusions.get(chrPair);
        List<FusionCohortData> fusionsByPosition = null;

        int fusionStartPos = fusion.JunctionPositions[SE_START];

        if(chrPairFusions == null)
        {
            chrPairFusions = Maps.newHashMap();
            mFusions.put(chrPair, chrPairFusions);

            fusionsByPosition = Lists.newArrayList();
            chrPairFusions.put(fusionStartPos, fusionsByPosition);
        }
        else
        {
            fusionsByPosition = chrPairFusions.get(fusionStartPos);

            if(fusionsByPosition == null)
            {
                fusionsByPosition = Lists.newArrayList();
                chrPairFusions.put(fusionStartPos, fusionsByPosition);
            }
            else
            {
                // check for a match
                FusionCohortData existingFusion = fusionsByPosition.stream().filter(x -> x.matches(fusion)).findFirst().orElse(null);

                if(existingFusion != null)
                {
                    existingFusion.addSample(sampleId, fusion.SplitFrags);
                    return;
                }
            }
        }

        FusionCohortData fusionCohortData = FusionCohortData.from(fusion);
        fusionCohortData.addSample(sampleId, fusion.SplitFrags);
        fusionsByPosition.add(fusionCohortData);
        ++mFusionCount;
    }

    private void writeFilteredFusion(final String sampleId, final List<FusionData> sampleFusions)
    {
        // mark passing fusions, and then include any which are related to them
        final List<FusionData> passingFusions = Lists.newArrayList();
        final List<FusionData> nonPassingFusionsWithRelated = Lists.newArrayList();

        for (FusionData fusion : sampleFusions)
        {
            if(fusion.Id == 3003911)
            {
                ISF_LOGGER.debug("here");
            }

            mFilters.markKnownGeneTypes(fusion);

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

        ISF_LOGGER.info("sample({}) writing {} passing fusions from total({} relatedToPass={})",
                sampleId, passingFusions.size(), sampleFusions.size(), relatedToPassing);

        try
        {
            final String outputFileName = mConfig.OutputDir + sampleId + ".isf." + PASS_FUSION_FILE_ID;
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(mFilteredFusionHeader);
            writer.write(",GeneTypeUp,GeneTypeDown");
            writer.newLine();

            for (FusionData fusion : passingFusions)
            {
                writer.write(fusion.rawData());
                writer.write(String.format(",%s,%s", fusion.getGeneTypes()[SE_START], fusion.getGeneTypes()[SE_END]));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion cohort file: {}", e.toString());
        }
    }

}
