package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.ALLELE_FREQUENCY;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.ANCHOR_DISTANCE;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.COHORT;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.FRAGMENT_COUNT;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_PAIR;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_PROM3;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.OTHER_PROM3;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_KNOWN;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_PROM3;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.hasKnownPairGene;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

import org.apache.commons.cli.CommandLine;

public class FusionFilters
{
    private final FusionCohortConfig mConfig;

    private final Map<String, Map<Integer, List<FusionCohortData>>> mCohortFusions;
    private final KnownFusionCache mKnownFusionCache;

    public FusionFilters(final FusionCohortConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mCohortFusions = Maps.newHashMap();

        if(mConfig.CohortFile != null)
        {
            loadCohortFile();
        }

        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(cmd);
    }

    private static final int KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS = 2;
    private static final int KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS = 4;
    private static final int MIN_ANCHOR_DISTANCE = 20;
    private static final double AF_KNOWN_SPLICED_TIER = 0.002;
    private static final double AF_KNOWN_TIER = 0.005;
    private static final double AF_UNKNOWN_TIER = 0.05;
    private static final int LOCAL_FUSION_THRESHOLD = 1000000;

    public static boolean isShortLocalFusion(final FusionData fusion)
    {
        return fusion.Chromosomes[SE_START].equals(fusion.Chromosomes[SE_END]) &&
                abs(fusion.JunctionPositions[SE_END] - fusion.JunctionPositions[SE_START]) < LOCAL_FUSION_THRESHOLD;
    }

    public boolean isPassingFusion(final FusionData fusion)
    {
        return isPassingFusion(fusion, fusion.getKnownFusionType(), fusion.hasKnownSpliceSites());
    }

    public boolean isPassingFusion(final FusionData fusion, final KnownGeneType knownType, boolean hasKnownSpliceSites)
    {
        if(knownType == KNOWN_PAIR && !isShortLocalFusion(fusion))
        {
            int requiredFragments = (hasKnownSpliceSites || hasKnownSpliceSite(fusion)) ?
                    KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS : KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS;

            if(fusion.totalFragments() < requiredFragments)
            {
                fusion.setFilter(FRAGMENT_COUNT);
                return false;
            }

            return true;
        }

        if(min(fusion.AnchorDistance[SE_START], fusion.AnchorDistance[SE_END]) < MIN_ANCHOR_DISTANCE && fusion.DiscordantFrags == 0)
        {
            fusion.setFilter(ANCHOR_DISTANCE);
            return false;
        }

        double requiredAF;
        int requiredFragments;

        if(hasKnownSpliceSites)
        {
            requiredAF = AF_KNOWN_SPLICED_TIER;
            requiredFragments = 2;
        }
        else if((fusion.JunctionTypes[SE_START] == FusionJunctionType.CANONICAL && fusion.JunctionTypes[SE_END] == FusionJunctionType.KNOWN)
            || (fusion.JunctionTypes[SE_START] == FusionJunctionType.KNOWN && fusion.JunctionTypes[SE_END] == FusionJunctionType.CANONICAL))
        {
            requiredAF = AF_KNOWN_TIER;
            requiredFragments = 3;
        }
        else if(fusion.JunctionTypes[SE_START] == FusionJunctionType.CANONICAL && fusion.JunctionTypes[SE_END] == FusionJunctionType.CANONICAL)
        {
            requiredAF = AF_KNOWN_TIER;
            requiredFragments = 4;
        }
        else
        {
            requiredAF = AF_UNKNOWN_TIER;
            requiredFragments = 10;
        }

        if(fusion.alleleFrequency() < requiredAF)
        {
            fusion.setFilter(ALLELE_FREQUENCY);
            return false;
        }

        if(fusion.totalFragments() < requiredFragments)
        {
            fusion.setFilter(FRAGMENT_COUNT);
            return false;
        }

        int cohortFreqLimit = hasKnownPairGene(knownType) ? 5 : 2;

        if(fusion.cohortFrequency() >= cohortFreqLimit)
        {
            fusion.setFilter(COHORT);
            return false;
        }

        return true;
    }

    private static boolean hasKnownSpliceSite(final FusionData fusion)
    {
        return fusion.hasKnownSpliceSites()
            || (fusion.JunctionTypes[SE_START] == FusionJunctionType.CANONICAL && fusion.JunctionTypes[SE_END] == FusionJunctionType.KNOWN)
            || (fusion.JunctionTypes[SE_START] == FusionJunctionType.KNOWN && fusion.JunctionTypes[SE_END] == FusionJunctionType.CANONICAL);
    }

    public static boolean hasSufficientKnownFusionFragments(final List<FusionData> fusions)
    {
        int knownSpliceSiteFragments = fusions.stream()
                .filter(x -> hasKnownSpliceSite(x))
                .mapToInt(x -> x.totalFragments()).sum();

        int totalFragments = fusions.stream().mapToInt(x -> x.totalFragments()).sum();

        return knownSpliceSiteFragments >= KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS || totalFragments >= KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS;
    }

    public FusionCohortData findCohortFusion(final FusionData fusion)
    {
        final String chrPair = formChromosomePair(fusion.Chromosomes[SE_START], fusion.Chromosomes[SE_END]);

        final Map<Integer,List<FusionCohortData>> chrPairFusions = mCohortFusions.get(chrPair);

        if(chrPairFusions == null)
            return null;

        final List<FusionCohortData> fusionsByPosition = chrPairFusions.get(fusion.JunctionPositions[SE_START]);

        if(fusionsByPosition == null)
            return null;

        return fusionsByPosition.stream().filter(x -> x.matches(fusion)).findFirst().orElse(null);
    }

    public void markKnownGeneTypes(final FusionData fusion)
    {
        boolean[] isKnown = {false, false};

        for(KnownFusionData knownFusionData : mKnownFusionCache.getDataByType(KnownFusionType.KNOWN_PAIR))
        {
            if(knownFusionData.FiveGene.equals(fusion.GeneNames[SE_START]) && knownFusionData.ThreeGene.equals(fusion.GeneNames[SE_END]))
            {
                fusion.setKnownFusionType(KNOWN_PAIR);
                return;
            }

            if(knownFusionData.FiveGene.equals(fusion.GeneNames[SE_START]))
                isKnown[SE_START] = true;

            if(knownFusionData.ThreeGene.equals(fusion.GeneNames[SE_END]))
                isKnown[SE_END] = true;
        }

        boolean[] isProm = { false, false };

        if(!isKnown[SE_START] && mKnownFusionCache.hasPromiscuousFiveGene(fusion.GeneNames[SE_START]))
            isProm[SE_START] = true;

        if(!isKnown[SE_END] && mKnownFusionCache.hasPromiscuousThreeGene(fusion.GeneNames[SE_END]))
            isProm[SE_END] = true;

        if(isKnown[SE_START])
        {
            if(isProm[SE_END])
                fusion.setKnownFusionType(KNOWN_PROM3);
            else
                fusion.setKnownFusionType(KNOWN_OTHER);
        }
        else if(isKnown[SE_END])
        {
            if(isProm[SE_START])
                fusion.setKnownFusionType(PROM5_KNOWN);
            else
                fusion.setKnownFusionType(KNOWN_OTHER);
        }
        else if(isProm[SE_START])
        {
            if(isProm[SE_END])
                fusion.setKnownFusionType(PROM5_PROM3);
            else
                fusion.setKnownFusionType(PROM5_OTHER);
        }
        else if(isProm[SE_END])
        {
            fusion.setKnownFusionType(OTHER_PROM3);
        }
        else
        {
            fusion.setKnownFusionType(OTHER);
        }

    }

    private void loadCohortFile()
    {
        ISF_LOGGER.info("loading cohort fusion file", mConfig.CohortFile);

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(mConfig.CohortFile));

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);

            lines.remove(0);

            int fusionCount = 0;
            for(String fusionData : lines)
            {
                FusionCohortData fusion = FusionCohortData.fromCsv(fusionData, fieldsMap);

                ++fusionCount;

                final String chrPair = formChromosomePair(fusion.Chromosomes[SE_START], fusion.Chromosomes[SE_END]);

                Map<Integer,List<FusionCohortData>> chrPairFusions = mCohortFusions.get(chrPair);
                List<FusionCohortData> fusionsByPosition = null;

                int fusionStartPos = fusion.JunctionPositions[SE_START];

                if(chrPairFusions == null)
                {
                    chrPairFusions = Maps.newHashMap();
                    mCohortFusions.put(chrPair, chrPairFusions);

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
                }

                fusionsByPosition.add(fusion);
            }

            ISF_LOGGER.info("loaded {} cohort fusions", fusionCount);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion cohort file({}): {}", mConfig.CohortFile.toString(), e.toString());
        }
    }
}
