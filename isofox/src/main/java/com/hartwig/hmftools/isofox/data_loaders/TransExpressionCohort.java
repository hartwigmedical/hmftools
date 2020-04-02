package com.hartwig.hmftools.isofox.data_loaders;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.TRANSCRIPT;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FIT_ALLOCATION;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;


import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class TransExpressionCohort
{
    private final DataLoaderConfig mConfig;
    private final Map<String,Integer> mFieldsMap;

    private final Map<Integer,TransExpressionData> mTranscriptExpressionData;

    private BufferedWriter mTransDistributionWriter;

    private static final int DISTRIBUTION_SIZE = 101;

    public TransExpressionCohort(final DataLoaderConfig config)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();
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
            // rewriteSampleTranscripts(sampleId, transDataList);
            ISF_LOGGER.debug("{}: sample({}) loaded transcript data", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples transcript files", mConfig.SampleData.SampleIds.size());

        writeTranscriptTpmPercentiles();

        // write a report for any re-occurring alt SJ
        // writeTranscriptFrequencyDistribution();
        closeBufferedWriter(mTransDistributionWriter);
    }

    public static void calcPercentileValues(final List<Double> tpmValues, final double[] percentileValues)
    {
        int sampleCount = tpmValues.size();

        // populate the upper and lower bounds
        double percSlots = percentileValues.length;

        double samplesPerPercentile = sampleCount/percSlots;

        for(int i = 0; i < percentileValues.length; ++i)
        {
            double lowerIndex = i * samplesPerPercentile;
            double upperIndex = lowerIndex + samplesPerPercentile * 0.9999;

            int lowerBound = (int)floor(lowerIndex);
            int upperBound = (int)ceil(upperIndex) - 1;
            upperBound = min(upperBound, tpmValues.size());

            if(lowerBound == upperBound)
            {
                percentileValues[i] = tpmValues.get(lowerBound);
                continue;
            }

            double tpmTotal = 0;
            double sampleTotal = 0;

            for(int s = lowerBound; s <= upperBound; ++s)
            {
                double tpm = tpmValues.get(s);

                double fractionOfTpm;

                if(s == lowerBound)
                {
                    fractionOfTpm = 1 - (lowerIndex - lowerBound);
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else if(s == upperBound)
                {
                    fractionOfTpm = upperIndex - upperBound;
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else
                {
                    ++sampleTotal;
                    tpmTotal += tpm;
                }
            }

            percentileValues[i] = tpmTotal / sampleTotal;
        }
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
            for(final Map.Entry<Integer,TransExpressionData> entry : mTranscriptExpressionData.entrySet())
            {
                final TransExpressionData expData = entry.getValue();

                final double[] percentileValues = new double[DISTRIBUTION_SIZE + 1];

                calcPercentileValues(expData.TpmValues, percentileValues);

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

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            int geneIdIndex = mFieldsMap.get(FLD_GENE_ID);
            int transIdIndex = mFieldsMap.get(FLD_TRANS_ID);
            int fitAllocIndex = mFieldsMap.get(FLD_FIT_ALLOCATION);
            int tpmIndex = mFieldsMap.get(FLD_TPM);

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                int transId = Integer.parseInt(items[transIdIndex]);

                TransExpressionData transExpData = mTranscriptExpressionData.get(transId);

                if(transExpData == null)
                {
                    transExpData = TransExpressionData.fromCsv(items, mFieldsMap);
                    mTranscriptExpressionData.put(transId, transExpData);
                }

                double fitAllocation = Double.parseDouble(items[fitAllocIndex]);
                double tmp = Double.parseDouble(items[tpmIndex]);

                transExpData.addSampleData(sampleId, fitAllocation, tmp);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load transcript data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

}
