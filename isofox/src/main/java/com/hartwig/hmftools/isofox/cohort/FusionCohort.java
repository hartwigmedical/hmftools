package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.fromCsv;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
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

    private int mFusionCount;

    public FusionCohort(final CohortConfig config)
    {
        mConfig = config;
        mFusions = Maps.newHashMap();
        mFieldsMap = Maps.newHashMap();
        mFusionCount = 0;
    }

    public void processFusionFiles()
    {
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

            final List<FusionCohortData> fusions = loadFile(fusionFile, sampleId);

            ISF_LOGGER.info("{}: sample({}) loaded {} fusion records", i, sampleId, fusions.size());
            totalProcessed += fusions.size();

            fusions.forEach(x -> addFusion(x, sampleId));

            if(mFusionCount > nextLog)
            {
                nextLog += 100000;
                ISF_LOGGER.info("total fusion count({})", mFusionCount);
            }
        }

        ISF_LOGGER.info("loaded {} fusion records, total fusion count({})", totalProcessed, mFusionCount);
        writeReoccurringFusions();
    }

    private List<FusionCohortData> loadFile(final Path filename, final String sampleId)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            return lines.stream()
                    .map(x -> FusionCohortData.fromCsv(x, mFieldsMap, sampleId))
                    .collect(Collectors.toList());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion file({}): {}", filename.toString(), e.toString());
            return null;
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

    private void writeReoccurringFusions()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("fusion_cohort.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            for(int se = SE_START; se <= SE_END; ++se)
            {
                String prefix = se == SE_START ? "Up" : "Down";

                if(se == SE_END)
                    writer.write(",");

                writer.write(String.format("GeneId%s", prefix));

                writer.write(String.format(",GeneName%s", prefix));
                writer.write(String.format(",Chr%s", prefix));
                writer.write(String.format(",Pos%s", prefix));
                writer.write(String.format(",Orient%s", prefix));
                writer.write(String.format(",JuncType%s", prefix));
            }

            writer.write(",SVType,SampleCount,TotalFragments,MaxFragments,Samples");
            writer.newLine();

            for(Map<Integer,List<FusionCohortData>> chrPairLists : mFusions.values())
            {
                for (List<FusionCohortData> fusionLists : chrPairLists.values())
                {
                    for (FusionCohortData fusion : fusionLists)
                    {
                        if (fusion.sampleCount() < mConfig.FusionMinSampleThreshold)
                            continue;

                        for (int se = SE_START; se <= SE_END; ++se)
                        {
                            if (se == SE_END)
                                writer.write(",");

                            writer.write(String.format("%s,%s,%s,%d,%d,%s",
                                    fusion.GeneIds[se], fusion.GeneNames[se], fusion.Chromosomes[se],
                                    fusion.JunctionPositions[se], fusion.JunctionOrientations[se], fusion.JunctionTypes[se]));
                        }

                        writer.write(String.format(",%s,%d,%d,%d",
                                fusion.SvType, fusion.sampleCount(), fusion.fragmentCount(), fusion.maxFragmentCount()));

                        writer.write(String.format(",%s", appendStrList(fusion.sampleIds(), ';')));

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
}
