package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.FUSION_SOURCE_ARRIBA;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ISOFOX_ID;

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

    public void compareFusions(final String sampleId, final List<FusionData> fusions, final List<FusionData> unfilteredFusions)
    {
        if(mWriter == null)
            return;

        final Map<String,List<ExternalFusionData>> externalFusionData = loadExternalFusionFiles(sampleId);

        final List<FusionData> unmatchedFusions = Lists.newArrayList();
        final List<ExternalFusionData> matchedExtFusions = Lists.newArrayList();
        int matchedFusionCount = 0;

        for(FusionData fusion : fusions)
        {
            List<ExternalFusionData> extFusions = Lists.newArrayList();

            // factor in orientations being the wrong way around
            List<ExternalFusionData> extFusions1 = externalFusionData.get(formChromosomePair(fusion.Chromosomes));

            if(extFusions1 != null)
                extFusions.addAll(extFusions1);

            if(!fusion.Chromosomes[SE_START].equals(fusion.Chromosomes[SE_END]))
            {
                List<ExternalFusionData> extFusions2 = externalFusionData.get(formChromosomePair(fusion.Chromosomes[SE_END], fusion.Chromosomes[SE_START]));
                if(extFusions2 != null)
                    extFusions.addAll(extFusions2);
            }

            if(extFusions == null || extFusions.isEmpty())
            {
                unmatchedFusions.add(fusion);
                continue;
            }

            ExternalFusionData extFusion = extFusions.stream().filter(x -> x.matches(fusion)).findFirst().orElse(null);

            if(extFusion == null)
            {
                unmatchedFusions.add(fusion);
            }
            else
            {
                ++matchedFusionCount;
                matchedExtFusions.add(extFusion);

                String matchType = extFusion.IsFiltered ? "MATCH" : "MATCH_EXT_UNFILTERED";

                writeMatchData(sampleId, matchType, fusion.Id, fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations,
                        fusion.SvType, fusion.GeneNames, fusion.SplitFrags + fusion.RealignedFrags, fusion.DiscordantFrags,
                        extFusion.SplitFragments, extFusion.DiscordantFragments,
                        extFusion.IsFiltered ? "" : extFusion.otherDataStr());
            }
        }

        for(FusionData fusion : unmatchedFusions)
        {
            writeMatchData(sampleId, "ISOFOX_ONLY", fusion.Id, fusion.Chromosomes, fusion.JunctionPositions,
                    fusion.JunctionOrientations, fusion.SvType, fusion.GeneNames, fusion.SplitFrags + fusion.RealignedFrags,
                    fusion.DiscordantFrags,0, 0, fusionId(fusion.Id));
        }

        int extUnmatchedCount = 0;
        for(List<ExternalFusionData> extFusionLists : externalFusionData.values())
        {
            for(ExternalFusionData extFusion : extFusionLists)
            {
                if(!extFusion.IsFiltered || matchedExtFusions.contains(extFusion))
                    continue;

                // check against isofox unfiltered
                if(!unfilteredFusions.isEmpty())
                {
                    FusionData unfilteredFusion = unfilteredFusions.stream().filter(x -> extFusion.matches(x)).findFirst().orElse(null);
                    if(unfilteredFusion != null)
                    {
                        writeMatchData(sampleId, "MATCH_ISF_UNFILTERED", unfilteredFusion.Id,
                                extFusion.Chromosomes, extFusion.JunctionPositions, extFusion.JunctionOrientations,
                                extFusion.SvType, extFusion.GeneNames,
                                unfilteredFusion.SplitFrags + unfilteredFusion.RealignedFrags, unfilteredFusion.DiscordantFrags,
                                extFusion.SplitFragments, extFusion.DiscordantFragments, "");

                        continue;
                    }
                }

                ++extUnmatchedCount;
                writeMatchData(sampleId, "EXT_ONLY", -1, extFusion.Chromosomes, extFusion.JunctionPositions, extFusion.JunctionOrientations,
                        extFusion.SvType, extFusion.GeneNames, 0, 0, extFusion.SplitFragments, extFusion.DiscordantFragments, "");
            }
        }

        ISF_LOGGER.info("sample({}) matched({}) from isofox fusions({}) external({})",
                sampleId, matchedFusionCount, fusions.size(), extUnmatchedCount + matchedFusionCount);
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("ext_fusions_compare.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,MatchType,FusionId,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown"
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

    private void writeMatchData(
            final String sampleId, final String matchType, int fusionsId, final String[] chromosomes, final int[] junctionPositions,
            final byte[] junctionOrientations, final String svType, final String[] geneNames, int junctFrags, int discFrags,
            int extJunctFrags, int extDiscFrags, final String otherData)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s,%d,%d,%d,%d",
                    sampleId, matchType, fusionsId >= 0 ? fusionId(fusionsId) : "NONE",
                    chromosomes[SE_START], chromosomes[SE_END],
                    junctionPositions[SE_START], junctionPositions[SE_END], junctionOrientations[SE_START], junctionOrientations[SE_END]));

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

    private final Map<String, List<ExternalFusionData>> loadExternalFusionFiles(final String sampleId)
    {
        final Map<String, List<ExternalFusionData>> mappedFusions = Maps.newHashMap();

        if(mConfig.Fusions.ComparisonSources.contains(FUSION_SOURCE_ARRIBA))
        {
            final String arribaMainFile = mConfig.RootDataDir + sampleId + ".fusions.tsv";

            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(arribaMainFile));
                lines.remove(0);

                List<ExternalFusionData> sampleFusions = lines.stream()
                        .map(x -> ExternalFusionData.loadArribaFusion(x))
                        .filter(x -> x != null)
                        .collect(Collectors.toList());

                if(mConfig.Fusions.CompareUnfiltered)
                {
                    final String arribaDiscardedFile = mConfig.RootDataDir + sampleId + ".fusions.discarded.tsv";
                    lines.clear();
                    lines.addAll(Files.readAllLines(Paths.get(arribaDiscardedFile)));
                    lines.remove(0);

                    List<ExternalFusionData> discardedFusions = lines.stream()
                            .map(x -> ExternalFusionData.loadArribaFusion(x))
                            .filter(x -> x != null)
                            .collect(Collectors.toList());

                    discardedFusions.forEach(x -> x.IsFiltered = false);
                    sampleFusions.addAll(discardedFusions);
                }

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
