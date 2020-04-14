package com.hartwig.hmftools.isofox.data_loaders;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.calcPercentileValues;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.TRANSCRIPT;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionData.TPM_COUNT;
import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionData.TPM_VALUE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_EFFECTIVE_LENGTH;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FITTED_FRAGMENTS;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class TransExpressionCohort
{
    private final DataLoaderConfig mConfig;

    private final Map<Integer,TransExpressionData> mTranscriptExpressionData;

    private BufferedWriter mTransDistributionWriter;

    private static final int DISTRIBUTION_SIZE = 101;

    public TransExpressionCohort(final DataLoaderConfig config)
    {
        mConfig = config;
        mTranscriptExpressionData = Maps.newHashMap();

        mTransDistributionWriter = null;
    }

    public void processTranscripts()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, TRANSCRIPT, filenames))
            return;

        initialiseWriter();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path transcriptsFile = filenames.get(i);

            loadFile(sampleId, transcriptsFile);
            ISF_LOGGER.debug("{}: sample({}) loaded transcript data", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples transcript files", mConfig.SampleData.SampleIds.size());

        writeTranscriptTpmPercentiles();

        closeBufferedWriter(mTransDistributionWriter);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("transcript_distribution.csv");
            mTransDistributionWriter = createBufferedWriter(outputFileName, false);

            mTransDistributionWriter.write("GeneId,GeneName,TransName");

            double distributionSize = DISTRIBUTION_SIZE;

            for(int i = 0; i <= DISTRIBUTION_SIZE; ++i)
            {
                mTransDistributionWriter.write(String.format(",Pct_%.2f", i/distributionSize));
            }

            mTransDistributionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private void writeTranscriptTpmPercentiles()
    {
        try
        {
            int sampleCount = mConfig.SampleData.SampleIds.size();

            for(final Map.Entry<Integer,TransExpressionData> entry : mTranscriptExpressionData.entrySet())
            {
                final TransExpressionData expData = entry.getValue();

                final double[] percentileValues = new double[DISTRIBUTION_SIZE + 1];

                calcPercentileValues(convertDistribution(expData.TpmValues, sampleCount), percentileValues);

                // skip a transcript if its 100th percentile TPM is below the threshold
                if(percentileValues[percentileValues.length - 1] < mConfig.TpmLogThreshold)
                    continue;

                mTransDistributionWriter.write(String.format("%s,%s,%s",
                        expData.GeneId, expData.GeneName, expData.TransName));

                for(int i = 0; i <= DISTRIBUTION_SIZE; ++i)
                {
                    mTransDistributionWriter.write(String.format(",%6.3e", percentileValues[i]));
                }

                mTransDistributionWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private void loadFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int transIdIndex = fieldsMap.get(FLD_TRANS_ID);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            boolean roundValues = mConfig.SampleData.SampleIds.size() >= 100;

            final Map<Integer,double[]> transTpmData = Maps.newHashMap();

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                int transId = Integer.parseInt(items[transIdIndex]);
                final String geneId = items[geneIdIndex];

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                TransExpressionData transExpData = mTranscriptExpressionData.get(transId);

                if (transExpData == null)
                {
                    transExpData = TransExpressionData.fromCsv(items, fieldsMap);
                    mTranscriptExpressionData.put(transId, transExpData);
                }

                double transPerMill = Double.parseDouble(items[tpmIndex]);

                if(roundValues)
                {
                    transPerMill = roundTPM(transPerMill, mConfig.TpmRounding);
                }

                transExpData.addSampleData(sampleId, transPerMill);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load transcript data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    public static double roundTPM(double tpm, double roundingFactor)
    {
        double scale = round(log10(tpm));
        double tick = pow(10, scale - roundingFactor);
        return round(tpm/tick) * tick;
    }

    public static List<Double> convertDistribution(final List<double[]> tpmValueCounts, int sampleCount)
    {
        final List<Double> tpmValues = Lists.newArrayListWithExpectedSize(sampleCount);

        for(int i = 0; i < tpmValueCounts.size(); ++i)
        {
            final double tpmData[] = tpmValueCounts.get(i);
            for(int j = 0; j < tpmData[TPM_COUNT]; ++j)
            {
                tpmValues.add(tpmData[TPM_VALUE]);
            }
        }

        return tpmValues;
    }
}
