package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.loadExternalFusionFiles;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;

public class FusionCohort
{
    private final CohortConfig mConfig;
    private final FusionCohortConfig mFusionConfig;

    private final FusionFilters mFilters;
    private final Map<String,Integer> mFieldsMap;

    // map of chromosome-pair to start position to list of fusion junctions
    private final Map<String, Map<Integer,List<FusionCohortData>>> mFusions;

    private final Map<String,List<ExternalFusionData>> mExternalFusionData;

    private int mFusionCount;

    /* Routines:
        1. generate a cohort file from multiple sample fusion files from Isofox
        2. write a set of filtered/passing fusions for each sample fusion file loaded
        3. run comparison of fitered/passing fusions between Isofox and external fusion files
    */

    public FusionCohort(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mFusionConfig = new FusionCohortConfig(cmd);
        mFilters = new FusionFilters(mFusionConfig);
        mFusions = Maps.newHashMap();
        mFieldsMap = Maps.newHashMap();
        mExternalFusionData = Maps.newHashMap();
        mFusionCount = 0;
    }

    public void processFusionFiles()
    {
        if(!mFusionConfig.GenerateCohort && mFusionConfig.ComparisonSources.isEmpty())
        {
            ISF_LOGGER.warn("no fusion functions configured");
            return;
        }

        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, FUSION, filenames))
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

            if(mFusionConfig.WriteFilteredFusions)
            {
                sampleFusions.stream().filter(x -> mFilters.isPassingFusion(x)).forEach(x -> writeFilteredFusion(x));
                continue;
            }

            // add to the fusion cohort collection
            sampleFusions.forEach(x -> addToCohortCache(x, sampleId));

            if(!mFusionConfig.ComparisonSources.isEmpty())
            {
                loadExternalFusionFiles(sampleId, mConfig, mFusionConfig, mExternalFusionData);
                compareExternalFusions();
                mFusions.clear();
            }

            if(mFusionCount > nextLog)
            {
                nextLog += 100000;
                ISF_LOGGER.info("total fusion count({})", mFusionCount);
            }
        }

        if(mFusionConfig.GenerateCohort)
        {
            ISF_LOGGER.info("loaded {} fusion records, total fusion count({})", totalProcessed, mFusionCount);
            FusionCohortData.writeCohortFusions(mFusions, mConfig, mFusionConfig);
            return;
        }
    }

    private void writeFilteredFusion(final FusionData fusion)
    {

    }

    private void compareExternalFusions()
    {
        /*
        #ISFOX PASES

        View(merge(merge(passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":") ,breakpoint2=paste(ChrDown,PosDown,sep=":")),arriba %>% distinct(breakpoint1,breakpoint2,filters,result),by=c('breakpoint1','breakpoint2'),all.x=T),
             arriba %>% distinct(breakpoint1,breakpoint2,filters,result) %>% mutate(revOrientation=TRUE),by.x=c('breakpoint1','breakpoint2'),by.y=c('breakpoint2','breakpoint1'),all.x=T) %>% select(CoverageUp,CoverageDown,TotalFragments,everything()))

        #ARRIBA PASSES
        View(merge(merge(passedArriba,passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":")),by=c('breakpoint1','breakpoint2'),all.x=T),
                   passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":"),revOrientation=TRUE),by.x=c('breakpoint1','breakpoint2'),by.y=c('breakpoint2','breakpoint1'),all.x=T))# %>%
           filter(is.na(FusionId.x),is.na(FusionId.y),!grepl("read-through",type)))
         */
    }

    private List<FusionData> loadSampleFile(final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);
            return lines.stream().map(x -> FusionData.fromCsv(x, mFieldsMap)).collect(Collectors.toList());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion file({}): {}", filename.toString(), e.toString());
            return Lists.newArrayList();
        }
    }

    private void addToCohortCache(final FusionData fusion, final String sampleId)
    {
        if(fusion.SplitFrags < mFusionConfig.MinFragCount)
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
}
