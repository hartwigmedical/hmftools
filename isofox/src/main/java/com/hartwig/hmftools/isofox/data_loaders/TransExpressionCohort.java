package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.TRANSCRIPT;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoader.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;
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

public class TransExpressionCohort
{
    private final DataLoaderConfig mConfig;

    private final Map<String,Integer> mFieldsMap;

    private BufferedWriter mTranscriptWriter;

    public TransExpressionCohort(final DataLoaderConfig config)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();
        mTranscriptWriter = null;
    }

    public void processTranscripts()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, TRANSCRIPT, filenames))
            return;

        initialiseWrite();

        int totalProcessed = 0;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path transcriptsFile = filenames.get(i);

            final List<TransExpressionData> transcripts = loadFile(transcriptsFile);

            ISF_LOGGER.debug("{}: sample({}) loaded {} transcript records", i, sampleId, transcripts.size());
            totalProcessed += transcripts.size();

            rewriteSampleTranscripts(sampleId, transcripts);
        }

        ISF_LOGGER.info("loaded {} transcript records", totalProcessed);

        // write a report for any re-occurring alt SJ
        // writeTranscriptFrequencyDistribution();
        closeBufferedWriter(mTranscriptWriter);
    }

    private void initialiseWrite()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("transcript_cohort.csv");
            mTranscriptWriter = createBufferedWriter(outputFileName, false);

            mTranscriptWriter.write("SampleId,GeneId,GeneName,TransName,EffectiveLength,FitAllocation,TPM");
            mTranscriptWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private void rewriteSampleTranscripts(final String sampleId, final List<TransExpressionData> transcripts)
    {
        // calculate a TPM for each transcript
        double totalFragsPerKb = transcripts.stream().mapToDouble(x -> x.fragsPerKb()).sum();
        double tmpFactor = totalFragsPerKb / 1e6;

        try
        {
            for(TransExpressionData transData : transcripts)
            {
                if (!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(transData.GeneId))
                    continue;

                transData.setTransPerM(transData.fragsPerKb()/tmpFactor);

                if(transData.transPerM() < mConfig.TpmLogThreshold)
                    continue;

                mTranscriptWriter.write(String.format("%s,%s,%s,%s,%d,%.1f,%.6f",
                        sampleId, transData.GeneId, transData.GeneName, transData.TransName,
                        transData.EffectiveLength, transData.FitAllocation, transData.transPerM()));
                mTranscriptWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private List<TransExpressionData> loadFile(final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            return lines.stream()
                    .map(x -> TransExpressionData.fromCsv(x, mFieldsMap))
                    .collect(Collectors.toList());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load transcript data file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

}
