package com.hartwig.hmftools.isofox.data_loaders;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.TRANSCRIPT;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_EFFECTIVE_LENGTH;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FIT_ALLOCATION;
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

        if(mConfig.ConsolidateTpmData)
            rewriteTranscriptTpmData();

        closeBufferedWriter(mTransDistributionWriter);
    }

    public static void calcPercentileValues(final List<Double> values, final double[] percentileValues)
    {
        int sampleCount = values.size();

        // populate the upper and lower bounds
        double percSlots = percentileValues.length;

        double samplesPerPercentile = sampleCount/percSlots;

        for(int i = 0; i < percentileValues.length; ++i)
        {
            double lowerIndex = i * samplesPerPercentile;
            double upperIndex = lowerIndex + samplesPerPercentile * 0.9999;

            int lowerBound = (int)floor(lowerIndex);
            int upperBound = (int)ceil(upperIndex) - 1;
            upperBound = min(upperBound, values.size());

            if(lowerBound == upperBound)
            {
                percentileValues[i] = values.get(lowerBound);
                continue;
            }

            double tpmTotal = 0;
            double sampleTotal = 0;

            for(int s = lowerBound; s <= upperBound; ++s)
            {
                double tpm = values.get(s);

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

    private void rewriteTranscriptTpmData()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("transcript_consolidated.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            String header = new StringJoiner(DELIMITER)
                    .add(FLD_SAMPLE_ID)
                    .add(FLD_GENE_ID)
                    .add(FLD_GENE_NAME)
                    .add(FLD_TRANS_NAME)
                    .add(FLD_FIT_ALLOCATION)
                    .add(FLD_EFFECTIVE_LENGTH)
                    .add(FLD_TPM)
                    .toString();

            writer.write(header);

            writer.newLine();

            for(final Map.Entry<Integer,TransExpressionData> entry : mTranscriptExpressionData.entrySet())
            {
                final TransExpressionData expData = entry.getValue();

                for(int i = 0; i < expData.SampleIds.size(); ++i)
                {
                    writer.write(String.format("%s,%s,%s,%s",
                            expData.SampleIds.get(i), expData.GeneId, expData.GeneName, expData.TransName));

                    writer.write(String.format(",%.0f,%d,%6.3e",
                            expData.FitAllocations.get(i), expData.EffectiveLengths.get(i), expData.TpmValues.get(i)));

                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
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
            int fitAllocIndex = fieldsMap.get(FLD_FIT_ALLOCATION);
            int effectiveLengthIndex = fieldsMap.get(FLD_EFFECTIVE_LENGTH);
            int tpmIndex = fieldsMap.containsKey(FLD_TPM) ? fieldsMap.get(FLD_TPM) : -1;

            final Map<Integer,double[]> transTpmData = Maps.newHashMap();

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                int transId = Integer.parseInt(items[transIdIndex]);
                final String geneId = items[geneIdIndex];

                boolean excludeGene = !mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId);

                TransExpressionData transExpData = null;

                if(!excludeGene)
                {
                    transExpData = mTranscriptExpressionData.get(transId);

                    if (transExpData == null)
                    {
                        transExpData = TransExpressionData.fromCsv(items, fieldsMap);
                        mTranscriptExpressionData.put(transId, transExpData);
                    }
                }

                // need all genes to calculate a TPM, and will filter out later
                if(tpmIndex >= 0 && excludeGene)
                    continue;

                double fitAllocation = Double.parseDouble(items[fitAllocIndex]);
                int effectiveLength = Integer.parseInt(items[effectiveLengthIndex]);

                if(tpmIndex >= 0)
                {
                    double transPerMill = tpmIndex >= 0 ? Double.parseDouble(items[tpmIndex]) : 0;
                    transExpData.addSampleData(sampleId, fitAllocation, transPerMill, effectiveLength);
                }
                else
                {
                    double[] tpmData = new double[TPM_TPM+1];
                    tpmData[TPM_FIT] = fitAllocation;
                    tpmData[TPM_LENGTH] = effectiveLength;
                    transTpmData.put(transId, tpmData);
                }
            }

            if(!transTpmData.isEmpty())
            {
                calcAndAddTpmData(sampleId, transTpmData);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load transcript data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private static final int TPM_FIT = 0;
    private static final int TPM_LENGTH = 1;
    private static final int TPM_FPK = 2;
    private static final int TPM_TPM = 3;

    private void calcAndAddTpmData(final String sampleId, final Map<Integer,double[]> transTpmData)
    {
        double totalFragsPerKb = 0;

        for(final double[] tpmData : transTpmData.values())
        {
            tpmData[TPM_FPK] = tpmData[TPM_LENGTH] > 0 ? tpmData[TPM_FIT] / (tpmData[TPM_LENGTH] / 1000.0) : 0;
            totalFragsPerKb += tpmData[TPM_FPK];
        }

        double tpmFactor = totalFragsPerKb / 1e6;

        for(Map.Entry<Integer,double[]> entry : transTpmData.entrySet())
        {
            int transId = entry.getKey();
            final double[] tpmData = entry.getValue();

            tpmData[TPM_TPM] = tpmData[TPM_FPK] / tpmFactor;

            TransExpressionData transExpData = mTranscriptExpressionData.get(transId);

            if(transExpData == null) // filtered out earlier
                continue;

            transExpData.addSampleData(sampleId, tpmData[TPM_FIT], tpmData[TPM_TPM], (int)tpmData[TPM_LENGTH]);
        }
    }

}
