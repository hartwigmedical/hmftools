package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.FUSION_SOURCE_ARRIBA;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class ExternalFusionCompare
{
    private final CohortConfig mConfig;
    private final BufferedWriter mWriter;

    public ExternalFusionCompare(final CohortConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("ext_fusions_compare.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("MatchType,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown"
                    + ",SVType,GeneNameUp,GeneNameDown,JuncFrags,ExtJuncFrags,DiscFrags,ExtDiscFrags,OtherData");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion comparison file: {}", e.toString());
            return null;
        }
    }

    public void compareFusions(final String sampleId, final List<FusionData> fusions)
    {
        if(mWriter == null)
            return;

        final Map<String,List<ExternalFusionData>> externalFusionData = loadExternalFusionFiles(sampleId);

        final List<FusionData> unmatchedFusions = Lists.newArrayList();
        int matchedFusionCount = 0;

        for(FusionData fusion : fusions)
        {
            List<ExternalFusionData> extFusions = externalFusionData.get(formChromosomePair(fusion.Chromosomes));

            if(extFusions == null)
            {
                unmatchedFusions.add(fusion);
                continue;
            }

            ExternalFusionData externalFusion = extFusions.stream().filter(x -> x.matches(fusion)).findFirst().orElse(null);

            if(externalFusion == null)
            {
                unmatchedFusions.add(fusion);
            }
            else
            {
                ++matchedFusionCount;
                writeMatchData(sampleId, "MATCH", fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations,
                        fusion.SvType, fusion.GeneNames, fusion.SplitFrags + fusion.RealignedFrags, fusion.DiscordantFrags,
                        externalFusion.SplitFragments, externalFusion.DiscordantFragments, "");


                extFusions.remove(externalFusion);
            }
        }

        for(FusionData fusion : unmatchedFusions)
        {
            writeMatchData(sampleId, "ISOFOX_ONLY", fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations,
                    fusion.SvType, fusion.GeneNames, fusion.SplitFrags + fusion.RealignedFrags, fusion.DiscordantFrags,
                    0, 0, "");
        }

        int extUnmatchedCount = 0;
        for(List<ExternalFusionData> extFusionLists : externalFusionData.values())
        {
            for(ExternalFusionData fusion : extFusionLists)
            {
                ++extUnmatchedCount;
                writeMatchData(sampleId, "EXT_ONLY", fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations,
                        fusion.SvType, fusion.GeneNames, 0, 0, fusion.SplitFragments, fusion.DiscordantFragments, "");
            }
        }

        ISF_LOGGER.info("sample({}) matched({}) from isofox fusions({}) external({})",
                sampleId, matchedFusionCount, fusions.size(), extUnmatchedCount + matchedFusionCount);

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

    private void writeMatchData(
            final String sampleId, final String matchType, final String[] chromosomes, final int[] junctionPositions,
            final byte[] junctionOrientations, final String svType, final String[] geneNames, int junctFrags, int discFrags,
            int extJunctFrags, int extDiscFrags, final String otherData)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%d,%d,%d,%d",
                    matchType, chromosomes[SE_START], chromosomes[SE_END], junctionPositions[SE_START], junctionPositions[SE_END],
                    junctionOrientations[SE_START], junctionOrientations[SE_END]));

            mWriter.write(String.format(",%s,%s,%s,%d,%d,%d,%d,%s",
                    svType, geneNames[SE_START], geneNames[SE_END], junctFrags, extJunctFrags, discFrags, extDiscFrags, otherData));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion comparison file: {}", e.toString());
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

    private final Map<String, List<ExternalFusionData>>  loadExternalFusionFiles(final String sampleId)
    {
        final Map<String, List<ExternalFusionData>> mappedFusions = Maps.newHashMap();

        if(mConfig.Fusions.ComparisonSources.contains(FUSION_SOURCE_ARRIBA))
        {
            final String arribaMainFile = mConfig.RootDataDir + sampleId + ".fusions.tsv";

            // final String arribaDiscardedFile = config.RootDataDir + sampleId + ".fusions.discarded.tsv";

            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(arribaMainFile));
                lines.remove(0);

                /*
                final List<String> lines2 = Files.readAllLines(Paths.get(arribaDiscardedFile));
                lines2.remove(0);
                lines.addAll(lines2);
                */

                List<ExternalFusionData> sampleFusions = lines.stream()
                        .map(x -> ExternalFusionData.loadArribaFusion(x))
                        .filter(x -> x != null)
                        .collect(Collectors.toList());

                for(ExternalFusionData fusion : sampleFusions)
                {
                    final String chrPair = formChromosomePair(fusion.Chromosomes);

                    List<ExternalFusionData> chrPairFusions = mappedFusions.get(chrPair);
                    if(chrPairFusions == null)
                        mappedFusions.put(chrPair, Lists.newArrayList(fusion));
                    else
                        chrPairFusions.add(fusion);
                }
            }
            catch(IOException e)
            {
                ISF_LOGGER.error("failed to load arriba fusion file({}): {}", arribaMainFile, e.toString());
            }
        }

        return mappedFusions;
    }
}
