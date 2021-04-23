package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class FusionCollection
{
    private final CohortConfig mConfig;

    // map of chromosome-pair to start position to list of fusion junctions
    private final Map<String, Map<Integer,List<FusionCohortData>>> mFusions;

    private int mFusionCount;
    private int mFusionsProcessed;
    private int mNextLog;

    private static final int LOG_FUSION_INTERVAL = 100000;

    public FusionCollection(final CohortConfig config)
    {
        mConfig = config;
        mFusions = Maps.newHashMap();
        mFusionCount = 0;
        mFusionsProcessed = 0;
        mNextLog = LOG_FUSION_INTERVAL;
    }

    public synchronized void addToCohortCache(final List<FusionData> fusions, final String sampleId)
    {
        fusions.forEach(x -> addToCohortCache(x, sampleId));
        mFusionsProcessed += fusions.size();

        if(mFusionCount > mNextLog)
        {
            mNextLog += LOG_FUSION_INTERVAL;
            ISF_LOGGER.info("total fusion count({})", mFusionCount);
        }
    }

    private void addToCohortCache(final FusionData fusion, final String sampleId)
    {
        if(fusion.supportingFragments() < mConfig.Fusions.MinFragCount)
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
                    existingFusion.addSample(sampleId, fusion.supportingFragments());
                    return;
                }
            }
        }

        FusionCohortData fusionCohortData = FusionCohortData.from(fusion);
        fusionCohortData.addSample(sampleId, fusion.supportingFragments());
        fusionsByPosition.add(fusionCohortData);
        ++mFusionCount;
    }

    public void writeCohortFusions()
    {
        ISF_LOGGER.info("loaded {} fusion records, total fusion count({})", mFusionsProcessed, mFusionCount);

        try
        {
            final String outputFileName = mConfig.formCohortFilename("fusion_cohort.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(FusionCohortData.header());
            writer.newLine();

            int cohortFusionCount = 0;

            for(Map<Integer, List<FusionCohortData>> chrPairLists : mFusions.values())
            {
                for (List<FusionCohortData> fusionLists : chrPairLists.values())
                {
                    for (FusionCohortData fusion : fusionLists)
                    {
                        if (fusion.sampleCount() < mConfig.Fusions.MinSampleThreshold)
                            continue;

                        ++cohortFusionCount;

                        writer.write(FusionCohortData.toCsv(fusion));
                        writer.newLine();
                    }
                }
            }

            ISF_LOGGER.info("wrote {} cohort fusion", cohortFusionCount);

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusion cohort file: {}", e.toString());
        }

    }

}
