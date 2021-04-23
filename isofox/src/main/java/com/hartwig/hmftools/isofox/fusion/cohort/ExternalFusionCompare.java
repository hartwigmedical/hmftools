package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilename;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.COL_BREAKEND_1;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.COL_BREAKEND_2;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.FUSION_SOURCE_ARRIBA;
import static com.hartwig.hmftools.isofox.fusion.cohort.ExternalFusionData.junctionMatch;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

public class ExternalFusionCompare
{
    private final CohortConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final BufferedWriter mWriter;
    private final List<String> mSampleResults;
    private final Map<String,Integer> mFieldsMap;

    private static final String MT_MATCH = "MATCH";
    private static final String MT_FILT_IN_ISOFOX = "FILTERED_IN_ISOFOX";
    private static final String MT_FILT_IN_EXT = "FILTERED_IN_";
    private static final String MT_ISOFOX_ONLY = "ISOFOX_ONLY";
    private static final String MT_EXT_ONLY = "_ONLY";

    public ExternalFusionCompare(final CohortConfig config, final BufferedWriter writer)
    {
        mConfig = config;

        mGeneTransCache = new EnsemblDataCache(config.EnsemblDataCache, RefGenomeVersion.V37);
        mGeneTransCache.setRequiredData(false, false, false, false);
        mGeneTransCache.setRequireGeneSynonyms();
        mGeneTransCache.load(true); // only need  genes, not transcript data

        mWriter = writer;
        mSampleResults = Lists.newArrayList();
        mFieldsMap = Maps.newHashMap();
    }

    public final List<String> getResults() { return mSampleResults; }

    public void compareFusions(final String sampleId, final List<FusionData> fusions)
    {
        if(mWriter == null)
            return;

        final Map<String,List<ExternalFusionData>> externalFusionData = loadExternalFusionData(sampleId);

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

                cacheMatchResults(sampleId, getMatchType(MT_MATCH), fusion.Id, fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations,
                        fusion.JunctionTypes, fusion.SvType, fusion.GeneIds, fusion.GeneNames, fusion.SplitFrags + fusion.RealignedFrags,
                        fusion.DiscordantFrags, extFusion.SplitFragments, extFusion.DiscordantFragments, fusion.CohortCount, "");
            }
        }

        if(!unmatchedFusions.isEmpty())
        {
            checkExternalUnfilteredFusionData(sampleId, unmatchedFusions);
        }

        final List<ExternalFusionData> unmatchedExtFusions = Lists.newArrayList();

        for(List<ExternalFusionData> extFusionList : externalFusionData.values())
        {
            for(ExternalFusionData extFusion : extFusionList)
            {
                // Arriba can have duplicates so de-deup for these
                if(matchedExtFusions.stream().anyMatch(x -> junctionMatch(
                        x.Chromosomes, extFusion.Chromosomes, x.JunctionPositions, extFusion.JunctionPositions))
                || unmatchedExtFusions.stream().anyMatch(x -> junctionMatch(
                        x.Chromosomes, extFusion.Chromosomes, x.JunctionPositions, extFusion.JunctionPositions)))
                {
                    continue;
                }

                unmatchedExtFusions.add(extFusion);
            }
        }

        if(!unmatchedExtFusions.isEmpty())
        {
            checkUnfilteredFusionData(sampleId, unmatchedExtFusions);
        }

        // writeResults();

        ISF_LOGGER.debug("sample({}) matched({}) from isofox fusions({}) external({})",
                sampleId, matchedFusionCount, fusions.size(), matchedFusionCount);
    }

    public static BufferedWriter initialiseWriter(final CohortConfig config)
    {
        try
        {
            final String outputFileName = config.formCohortFilename("ext_fusions_compare.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,MatchType,FusionId,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown,JuncTypeUp,JuncTypeDown"
                    + ",SVType,GeneIdUp,GeneIdDown,GeneNameUp,GeneNameDown,JuncFrags,ExtJuncFrags,DiscFrags,ExtDiscFrags,CohortCount,OtherData");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion comparison file: {}", e.toString());
            return null;
        }
    }

    private void cacheMatchResults(
            final String sampleId, final String matchType, int fusionsId, final String[] chromosomes, final int[] junctionPositions,
            final byte[] junctionOrientations, final FusionJunctionType[] junctionTypes, final String svType, final String[] geneIds,
            final String[] geneNames, int junctFrags, int discFrags, int extJunctFrags, int extDiscFrags, int cohortCount, final String otherData)
    {
        String matchResultStr = String.format("%s,%s,%s,%s,%s,%d,%d,%d,%d,%s,%s",
                sampleId, matchType, fusionsId >= 0 ? fusionId(fusionsId) : "NONE",
                chromosomes[SE_START], chromosomes[SE_END], junctionPositions[SE_START], junctionPositions[SE_END],
                junctionOrientations[SE_START], junctionOrientations[SE_END], junctionTypes[SE_START], junctionTypes[SE_END]);

        matchResultStr += String.format(",%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%s",
                svType, geneIds[SE_START], geneIds[SE_END], geneNames[SE_START], geneNames[SE_END],
                junctFrags, extJunctFrags, discFrags, extDiscFrags, cohortCount, otherData);

        mSampleResults.add(matchResultStr);
    }

    public static void writeResults(final BufferedWriter writer, final List<String> results)
    {
        try
        {
            for(String result : results)
            {
                writer.write(result);
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion comparison file: {}", e.toString());
        }
    }

    private void checkUnfilteredFusionData(
            final String sampleId, final List<ExternalFusionData> unmatchedExternalFusions)
    {
        if(mConfig.Fusions.CompareUnfiltered)
        {
            final String unfilteredFile = formSampleFilename(mConfig, sampleId, FUSION);

            try
            {
                BufferedReader fileReader = new BufferedReader(new FileReader(unfilteredFile));

                String line = fileReader.readLine(); // skip header

                if(mFieldsMap.isEmpty())
                {
                    mFieldsMap.putAll(createFieldsIndexMap(line, DELIMITER));
                }

                while ((line = fileReader.readLine()) != null)
                {
                    // extract the required data from the unfiltered fusion record if a match is found
                    final FusionData unfilteredFusion = FusionData.fromCsv(line, mFieldsMap);

                    int index = 0;
                    while(index < unmatchedExternalFusions.size())
                    {
                        ExternalFusionData extFusion = unmatchedExternalFusions.get(index);

                        if(junctionMatch(
                                extFusion.Chromosomes, unfilteredFusion.Chromosomes, extFusion.JunctionPositions, unfilteredFusion.JunctionPositions))
                        {
                            final String otherData = String.format("reason=%s anchorDistance=%d af=%.3f",
                                    unfilteredFusion.getFilter().toString(),
                                    min(unfilteredFusion.AnchorDistance[SE_START], unfilteredFusion.AnchorDistance[SE_END]),
                                    unfilteredFusion.alleleFrequency());

                            cacheMatchResults(sampleId, getMatchType(MT_FILT_IN_ISOFOX), unfilteredFusion.Id,
                                    extFusion.Chromosomes, extFusion.JunctionPositions, extFusion.JunctionOrientations,
                                    unfilteredFusion.JunctionTypes, extFusion.SvType, unfilteredFusion.GeneIds, extFusion.GeneNames,
                                    unfilteredFusion.supportingFragments(), unfilteredFusion.DiscordantFrags,
                                    extFusion.SplitFragments, extFusion.DiscordantFragments, unfilteredFusion.CohortCount, otherData);

                            unmatchedExternalFusions.remove(index);

                            if(unmatchedExternalFusions.isEmpty())
                                return;

                            continue;
                        }

                        ++index;
                    }
                }
            }
            catch(IOException e)
            {
                ISF_LOGGER.error("failed to load unfiltered fusion file({}): {}", unfilteredFile, e.toString());
            }
        }

        for(ExternalFusionData extFusion : unmatchedExternalFusions)
        {
            // look up geneId data in Ensembl cache, including in synonymns if not found
            final String[] geneIds = {"", ""};
            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(extFusion.GeneNames[fs]);

                if(geneData == null)
                {
                    geneData = mGeneTransCache.getGeneDataBySynonym(extFusion.GeneNames[fs], extFusion.Chromosomes[fs]);
                }

                if(geneData != null)
                {
                    geneIds[fs] = geneData.GeneId;
                }
            }

            cacheMatchResults(sampleId, getMatchType(MT_EXT_ONLY), -1, extFusion.Chromosomes, extFusion.JunctionPositions,
                extFusion.JunctionOrientations, extFusion.JunctionTypes, extFusion.SvType, geneIds, extFusion.GeneNames,
                0, 0, extFusion.SplitFragments, extFusion.DiscordantFragments, 0, "");
        }

    }

    private String getMatchType(final String matchStr)
    {
        if(matchStr.equals(MT_MATCH) || matchStr.equals(MT_ISOFOX_ONLY) || matchStr.equals(MT_FILT_IN_ISOFOX))
            return matchStr;

        if(matchStr.equals(MT_EXT_ONLY))
            return mConfig.Fusions.ComparisonSource + matchStr;
        else
            return matchStr + mConfig.Fusions.ComparisonSource;
    }

    private Map<String, List<ExternalFusionData>> loadExternalFusionData(final String sampleId)
    {
        final Map<String, List<ExternalFusionData>> mappedFusions = Maps.newHashMap();

        if(mConfig.Fusions.ComparisonSource.equals(FUSION_SOURCE_ARRIBA))
        {
            final String arribaFile = mConfig.RootDataDir + sampleId + ".fusions.tsv";

            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(arribaFile));
                lines.remove(0);

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
                ISF_LOGGER.error("failed to load arriba fusion file({}): {}", arribaFile, e.toString());
            }
        }

        return mappedFusions;
    }

    private void checkExternalUnfilteredFusionData(final String sampleId, final List<FusionData> fusions)
    {
        if(mConfig.Fusions.CompareUnfiltered && mConfig.Fusions.ComparisonSource.equals(FUSION_SOURCE_ARRIBA))
        {
            final String arribaDiscardedFile = mConfig.RootDataDir + sampleId + ".fusions.discarded.tsv";

            try
            {
                BufferedReader fileReader = new BufferedReader(new FileReader(arribaDiscardedFile));

                String line = fileReader.readLine(); // skip header

                while ((line = fileReader.readLine()) != null)
                {
                    // extract the required data from the unfiltered fusion record if a match is found
                    final String[] items = line.split("\t", -1);

                    final String[] chromosomes = { items[COL_BREAKEND_1].split(":")[0], items[COL_BREAKEND_2].split(":")[0] };

                    if(!HumanChromosome.contains(chromosomes[SE_START]) || !HumanChromosome.contains(chromosomes[SE_END]))
                        continue;

                    final int[] positions = {
                            Integer.parseInt(items[COL_BREAKEND_1].split(":")[1]),
                            Integer.parseInt(items[COL_BREAKEND_2].split(":")[1]) };

                    int index = 0;
                    while(index < fusions.size())
                    {
                        FusionData fusion = fusions.get(index);

                        if(junctionMatch(fusion.Chromosomes, chromosomes, fusion.JunctionPositions, positions))
                        {
                            int splitFragments = Integer.parseInt(items[11]) + Integer.parseInt(items[12]);
                            int discordantFragments = Integer.parseInt(items[13]);
                            final String otherData = String.format("Conf=%s;Filters=%s",
                                    items[16], items[19].replaceAll(",",";"));

                            cacheMatchResults(sampleId, getMatchType(MT_FILT_IN_EXT), fusion.Id, fusion.Chromosomes,
                                    fusion.JunctionPositions, fusion.JunctionOrientations, fusion.JunctionTypes, fusion.SvType,
                                    fusion.GeneIds, fusion.GeneNames, fusion.supportingFragments(), fusion.DiscordantFrags,
                                    splitFragments, discordantFragments, fusion.CohortCount, otherData);

                            fusions.remove(index);

                            if(fusions.isEmpty())
                                return;

                            continue;
                        }

                        ++index;
                    }
                }
            }
            catch(IOException e)
            {
                ISF_LOGGER.error("failed to load Arriba discarded fusion file({}): {}", arribaDiscardedFile, e.toString());
            }
        }

        for(FusionData fusion : fusions)
        {
            cacheMatchResults(sampleId, getMatchType(MT_ISOFOX_ONLY), fusion.Id, fusion.Chromosomes, fusion.JunctionPositions,
                    fusion.JunctionOrientations, fusion.JunctionTypes, fusion.SvType, fusion.GeneIds, fusion.GeneNames,
                    fusion.SplitFrags + fusion.RealignedFrags, fusion.DiscordantFrags, 0, 0,
                    fusion.CohortCount, fusionId(fusion.Id));
        }
    }

}
