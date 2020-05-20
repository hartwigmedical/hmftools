package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.KnownFusionData.FIVE_GENE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.THREE_GENE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionData.FILTER_COHORT;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionData.FILTER_SUPPORT;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_PAIR;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.KNOWN_PROM3;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.OTHER_PROM3;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_KNOWN;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_OTHER;
import static com.hartwig.hmftools.isofox.fusion.cohort.KnownGeneType.PROM5_PROM3;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

import org.apache.commons.cli.CommandLine;

public class FusionFilters
{
    private final FusionCohortConfig mConfig;

    private final Map<String, Map<Integer, List<FusionCohortData>>> mCohortFusions;
    private final KnownFusionData mKnownFusionData;

    public FusionFilters(final FusionCohortConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mCohortFusions = Maps.newHashMap();

        if(mConfig.CohortFile != null)
        {
            loadCohortFile();
        }

        mKnownFusionData = new KnownFusionData();
        mKnownFusionData.loadFromFile(cmd);
    }

    private static final int MIN_ANCHOR_DISTANCE = 20;
    private static final double AF_KNOWN_TIER = 0.005;
    private static final double AF_UNKNOWN_TIER = 0.05;

    public boolean isPassingFusion(final FusionData fusion)
    {

        /*
        - AlleleFrequency = SplitFrags+RealignedFrags)/pmax(CoverageUp,CoverageDown)
        - Length = if not BND then abs(PosUp-PosDown)
        - not in PON until in a Known/Promiscuous category above
        - MaxAnchorLengthUp>20&MaxAnchorLengthDown>20)|DiscordantFrags>0),
        - (JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN'&AF>=0.005&TotalFragments>=2)|  #splice site - splice site
      (((JuncTypeUp=='CANONICAL'&JuncTypeDown=='KNOWN')|(JuncTypeUp=='KNOWN'&JuncTypeDown=='CANONICAL'))&AF>=0.005&TotalFragments>=3)|  #canonical - splice
      (JuncTypeUp=='CANONICAL'&JuncTypeDown=='CANONICAL'&AF>=0.005&TotalFragments>=4)|  #canonical - canonical
      (AF>=0.05&TotalFragments>=10)) %>% select(len,SVType,everything()))
         */

        if(fusion.getKnownFusionType() == KNOWN_PAIR)
            return true;

        if(min(fusion.AnchorDistance[SE_START], fusion.AnchorDistance[SE_END]) < MIN_ANCHOR_DISTANCE && fusion.DiscordantFrags == 0)
        {
            fusion.setFilter(FILTER_SUPPORT);
            return false;
        }

        double requiredAF;
        int requiredFragments;

        if(fusion.hasKnownSpliceSites())
        {
            requiredAF = AF_KNOWN_TIER;
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
            fusion.setFilter(FILTER_SUPPORT);
            return false;
        }

        if(fusion.totalFragments() < requiredFragments)
        {
            fusion.setFilter(FILTER_SUPPORT);
            return false;
        }

        if(fusion.cohortFrequency() > 0 || matchesCohortFusion(fusion))
        {
            fusion.setFilter(FILTER_COHORT);
            return false;
        }

        return true;
    }

    public boolean matchesCohortFusion(final FusionData fusion)
    {
        return findCohortFusion(fusion) != null;
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

        for(final String[] knownPair : mKnownFusionData.knownPairs())
        {
            if(knownPair[FIVE_GENE].equals(fusion.GeneNames[SE_START]) && knownPair[THREE_GENE].equals(fusion.GeneNames[SE_END]))
            {
                fusion.setKnownFusionType(KNOWN_PAIR);
                return;
            }

            if(knownPair[FIVE_GENE].equals(fusion.GeneNames[SE_START]))
                isKnown[SE_START] = true;

            if(knownPair[THREE_GENE].equals(fusion.GeneNames[SE_END]))
                isKnown[SE_END] = true;
        }

        boolean[] isProm = { false, false };

        if(!isKnown[SE_START] && mKnownFusionData.promiscuousFiveGenes().contains(fusion.GeneNames[SE_START]))
            isProm[SE_START] = true;

        if(!isKnown[SE_END] && mKnownFusionData.promiscuousThreeGenes().contains(fusion.GeneNames[SE_END]))
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
