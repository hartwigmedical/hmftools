package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FusionFilters
{
    private final FusionCohortConfig mConfig;

    private final Map<String, Map<Integer, List<FusionCohortData>>> mCohortFusions;

    public FusionFilters(final FusionCohortConfig config)
    {
        mConfig = config;

        mCohortFusions = Maps.newHashMap();

        if(mConfig.CohortFile != null)
        {
            loadCohortFile();
        }

    }

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


        return true;
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
