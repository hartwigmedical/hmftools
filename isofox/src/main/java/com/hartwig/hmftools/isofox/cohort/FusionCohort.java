package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.cohort.ExternalFusionData.FUSION_SOURCE_ARRIBA;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FusionCohort
{
    private final CohortConfig mConfig;

    private final Map<String,Integer> mFieldsMap;

    // map of chromosome-pair to start position to list of fusion junctions
    private final Map<String, Map<Integer,List<FusionCohortData>>> mFusions;
    private final Map<String, Map<Integer,List<FusionCohortData>>> mCohortFusions;

    private final Map<String,List<ExternalFusionData>> mExternalFusionData;

    private int mFusionCount;

    public FusionCohort(final CohortConfig config)
    {
        mConfig = config;
        mFusions = Maps.newHashMap();
        mCohortFusions = Maps.newHashMap();
        mFieldsMap = Maps.newHashMap();
        mExternalFusionData = Maps.newHashMap();
        mFusionCount = 0;
    }

    public void processFusionFiles()
    {
        if(!mConfig.FusionGenerateCohort && mConfig.FusionComparisonSources.isEmpty())
        {
            ISF_LOGGER.warn("no fusion functions configured");
            return;
        }

        if(mConfig.FusionCohortFile != null)
        {
            ISF_LOGGER.info("loading cohort fusion file", mConfig.FusionCohortFile);
            loadCohortFile(mConfig.FusionCohortFile);
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

            int sampleFusions = loadSampleFile(fusionFile, sampleId);

            ISF_LOGGER.info("{}: sample({}) loaded {} fusions", i, sampleId, sampleFusions);

            if(!mConfig.FusionComparisonSources.isEmpty())
            {
                loadExternalFusionFiles(sampleId);
                compareExternalFusions();
                mFusions.clear();
            }

            if(mFusionCount > nextLog)
            {
                nextLog += 100000;
                ISF_LOGGER.info("total fusion count({})", mFusionCount);
            }
        }

        if(mConfig.FusionGenerateCohort)
        {
            ISF_LOGGER.info("loaded {} fusion records, total fusion count({})", totalProcessed, mFusionCount);
            writeCohortFusions();
            return;
        }
    }

    private void compareExternalFusions()
    {
        /*
        passFusions = (fusions %>% mutate(AF=(SplitFrags+RealignedFrags)/pmax(CoverageUp,CoverageDown),len=ifelse(SVType!='BND',abs(PosDown-PosUp),0)) %>% filter(
      #!(grepl('HLA',GeneNameUp)),!(grepl('HLA',GeneNameDown)), #HLA issues (presumably sorted out by PON)
      #!(grepl('IGH',GeneNameUp)),!(grepl('IGK',GeneNameUp)),!(grepl('IGH',GeneNameDown)),!(grepl('IGK',GeneNameDown)),!(grepl('IGL',GeneNameUp)),!(grepl('IGL',GeneNameDown)), #HLA issues (presumably sorted out by PON)
      ((MaxAnchorLengthUp>20&MaxAnchorLengthDown>20)|DiscordantFrags>0),
      #(len>500000|(JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN')|SVType=='INV'|SVType=='BND'),
      (JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN'&AF>=0.005&TotalFragments>=2)|  #splice site - splice site
      (((JuncTypeUp=='CANONICAL'&JuncTypeDown=='KNOWN')|(JuncTypeUp=='KNOWN'&JuncTypeDown=='CANONICAL'))&AF>=0.005&TotalFragments>=3)|  #canonical - splice
      (JuncTypeUp=='CANONICAL'&JuncTypeDown=='CANONICAL'&AF>=0.005&TotalFragments>=4)|  #canonical - canonical
      (AF>=0.05&TotalFragments>=10)) %>% select(len,SVType,everything()))
    passFusions=(merge(passFusions,cohortFusions %>% select(PosUp,OrientUp,PosDown,OrientDown,SampleCount),by=c('PosUp','OrientUp','PosDown','OrientDown'),all.x=T)) %>% filter((GeneNameUp=='TMPRSS2'|is.na(SampleCount))) ### TMPRSS2 filter temporary

    #ISFOX PASES

    View(merge(merge(passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":") ,breakpoint2=paste(ChrDown,PosDown,sep=":")),arriba %>% distinct(breakpoint1,breakpoint2,filters,result),by=c('breakpoint1','breakpoint2'),all.x=T),
         arriba %>% distinct(breakpoint1,breakpoint2,filters,result) %>% mutate(revOrientation=TRUE),by.x=c('breakpoint1','breakpoint2'),by.y=c('breakpoint2','breakpoint1'),all.x=T) %>% select(CoverageUp,CoverageDown,TotalFragments,everything()))

    #ARRIBA PASSES
    View(merge(merge(passedArriba,passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":")),by=c('breakpoint1','breakpoint2'),all.x=T),
               passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":"),revOrientation=TRUE),by.x=c('breakpoint1','breakpoint2'),by.y=c('breakpoint2','breakpoint1'),all.x=T))# %>%
       filter(is.na(FusionId.x),is.na(FusionId.y),!grepl("read-through",type)))
         */
    }

    private int loadSampleFile(final Path filename, final String sampleId)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            int fusionCount = 0;
            for(String fusionData : lines)
            {
                FusionCohortData fusion = FusionCohortData.fromSampleFusion(fusionData, mFieldsMap, sampleId);
                ++fusionCount;
                addFusion(fusion, sampleId);
            }

            return fusionCount;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion file({}): {}", filename.toString(), e.toString());
            return 0;
        }
    }

    private void addFusion(final FusionCohortData fusion, final String sampleId)
    {
        if(fusion.fragmentCount() < mConfig.FusionMinFragCount)
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

        Map<Integer,List<FusionCohortData>> fusionsByPosition = mFusions.get(chrPair);

        int fusionStartPos = fusion.JunctionPositions[SE_START];

        if(fusionsByPosition == null)
        {
            fusionsByPosition = Maps.newHashMap();
            fusionsByPosition.put(fusionStartPos, Lists.newArrayList(fusion));
            mFusions.put(chrPair, fusionsByPosition);
            ++mFusionCount;
            return;
        }

        List<FusionCohortData> fusions = fusionsByPosition.get(fusionStartPos);

        if(fusions == null)
        {
            fusionsByPosition.put(fusionStartPos, Lists.newArrayList(fusion));
            ++mFusionCount;
            return;
        }

        // check for a match
        FusionCohortData existingFusion = fusions.stream().filter(x -> x.matches(fusion)).findFirst().orElse(null);

        if(existingFusion != null)
        {
            existingFusion.addSample(sampleId, fusion.fragmentCount());
        }
        else
        {
            ++mFusionCount;
            fusions.add(fusion);
        }
    }

    private void writeCohortFusions()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("fusion_cohort.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(FusionCohortData.cohortFusionHeader());
            writer.newLine();

            for(Map<Integer,List<FusionCohortData>> chrPairLists : mFusions.values())
            {
                for (List<FusionCohortData> fusionLists : chrPairLists.values())
                {
                    for (FusionCohortData fusion : fusionLists)
                    {
                        if (fusion.sampleCount() < mConfig.FusionMinSampleThreshold)
                            continue;

                        writer.write(fusion.toCohortFusion());
                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion cohort file: {}", e.toString());
        }
    }

    private void loadExternalFusionFiles(final String sampleId)
    {
        mExternalFusionData.clear();

        if(mConfig.FusionComparisonSources.contains(FUSION_SOURCE_ARRIBA))
        {
            final String arribaMainFile = mConfig.RootDataDir + sampleId + ".fusions.tsv";
            final String arribaDiscardedFile = mConfig.RootDataDir + sampleId + ".fusions.discarded.tsv";

            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(arribaMainFile));
                lines.remove(0);

                final List<String> lines2 = Files.readAllLines(Paths.get(arribaDiscardedFile));
                lines2.remove(0);
                lines.addAll(lines2);

                List<ExternalFusionData> fusions = lines.stream()
                        .map(x -> ExternalFusionData.loadArribaFusion(x))
                        .filter(x -> x != null)
                        .collect(Collectors.toList());

                for(ExternalFusionData fusion : fusions)
                {
                    final String chrPair = formChromosomePair(fusion.Chromosomes);

                    List<ExternalFusionData> chrPairFusions = mExternalFusionData.get(chrPair);
                    if(chrPairFusions == null)
                        mExternalFusionData.put(chrPair, Lists.newArrayList(fusion));
                    else
                        chrPairFusions.add(fusion);
                }
            }
            catch(IOException e)
            {
                ISF_LOGGER.error("failed to load arriba fusion file({}): {}", arribaMainFile, e.toString());
                return;
            }
        }
    }

    private void loadCohortFile(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);

            lines.remove(0);

            int fusionCount = 0;
            for(String fusionData : lines)
            {
                FusionCohortData fusion = FusionCohortData.fromCohortFusion(fusionData, fieldsMap);

                ++fusionCount;

                final String chrPair = formChromosomePair(fusion.Chromosomes[SE_START], fusion.Chromosomes[SE_END]);

                Map<Integer,List<FusionCohortData>> fusionsByPosition = mCohortFusions.get(chrPair);

                int fusionStartPos = fusion.JunctionPositions[SE_START];

                if(fusionsByPosition == null)
                {
                    fusionsByPosition = Maps.newHashMap();
                    fusionsByPosition.put(fusionStartPos, Lists.newArrayList(fusion));
                    mCohortFusions.put(chrPair, fusionsByPosition);
                    continue;
                }

                List<FusionCohortData> fusions = fusionsByPosition.get(fusionStartPos);

                if(fusions == null)
                    fusionsByPosition.put(fusionStartPos, Lists.newArrayList(fusion));
                else
                    fusions.add(fusion);
            }

            ISF_LOGGER.info("loaded {} cohort fusions", fusionCount);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion cohort file({}): {}", filename.toString(), e.toString());
        }
    }

}
