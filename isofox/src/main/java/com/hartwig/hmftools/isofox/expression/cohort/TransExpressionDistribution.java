package com.hartwig.hmftools.isofox.expression.cohort;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.TRANSCRIPT_DISTRIBUTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionData.RATE_COUNT;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionData.RATE_VALUE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class TransExpressionDistribution
{
    private final CohortConfig mConfig;

    private final Map<Integer,TransExpressionData> mTranscriptExpressionData;

    private BufferedWriter mTransDistributionWriter;

    private final Map<String,double[]> mCohortTranscriptDistribution;
    private final Map<String,double[]> mCancerTypeTranscriptDistribution;

    public static final int DISTRIBUTION_SIZE = PERCENTILE_COUNT; // percentiles from 0 to 100

    public TransExpressionDistribution(final CohortConfig config)
    {
        mConfig = config;
        mTranscriptExpressionData = Maps.newHashMap();
        mCohortTranscriptDistribution = Maps.newHashMap();
        mCancerTypeTranscriptDistribution = Maps.newHashMap();

        mTransDistributionWriter = null;

        if(mConfig.CohortTransFile != null)
        {
            loadCohortDistribution(
                    mConfig.CohortTransFile, mCohortTranscriptDistribution,
                    "transcript", DISTRIBUTION_SIZE + 1, Lists.newArrayList());
        }

        if(mConfig.CancerTransFile != null)
        {
            loadCohortDistribution(mConfig.CancerTransFile, mCancerTypeTranscriptDistribution,
                    "transcript", DISTRIBUTION_SIZE + 1, Lists.newArrayList());
        }
    }

    public void processSampleTranscriptFiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, TRANSCRIPT_DISTRIBUTION, filenames))
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

            for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
            {
                mTransDistributionWriter.write(String.format(",Pct_%d", i));
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

                final double[] percentileValues = new double[DISTRIBUTION_SIZE];

                final List<Double> tmpValues = convertDistribution(expData.TpmValues, sampleCount);
                final double[] tmpDataList = convertList(tmpValues);
                calcPercentileValues(tmpDataList, percentileValues);

                // skip a transcript if its 100th percentile TPM is below the threshold
                if(percentileValues[percentileValues.length - 1] < mConfig.TpmLogThreshold)
                    continue;

                mTransDistributionWriter.write(String.format("%s,%s,%s",
                        expData.GeneId, expData.GeneName, expData.TransName));

                for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
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
            for(int j = 0; j < tpmData[RATE_COUNT]; ++j)
            {
                tpmValues.add(tpmData[RATE_VALUE]);
            }
        }

        return tpmValues;
    }

    public static void loadCohortDistribution(
            final String inputFile, final Map<String,double[]> percentilesMap,
            final String fileType, int expectedColCount, final List<String> restrictions)
    {
        if(!Files.exists(Paths.get(inputFile)))
        {
            ISF_LOGGER.error("invalid cohort {} distribution file", fileType);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(inputFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty {} distribution file({})", fileType, inputFile);
                return;
            }

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                if (items.length != expectedColCount)
                {
                    ISF_LOGGER.error("invalid {} distribution data length({}) vs expected({}): {}",
                            fileType, items.length, expectedColCount, line);
                    return;
                }

                final String itemName = items[0];

                if(!restrictions.isEmpty() && !restrictions.contains(itemName))
                    continue;

                double[] percentileData = new double[DISTRIBUTION_SIZE];

                int startIndex = fileType.equals("gene") ? 2 : 1;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double tpm = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = tpm;
                }

                percentilesMap.put(itemName, percentileData);
            }

            ISF_LOGGER.info("loaded {} distribution from file({})", fileType, inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load {} distribution file({}): {}", fileType, inputFile, e.toString());
        }
    }

    public static double getTpmMedian(final Map<String,double[]> transPercentilesMap, final String transName)
    {
        final double[] transPercentiles = transPercentilesMap.get(transName);

        if(transPercentiles == null)
            return -1;

        return (transPercentiles[49] + transPercentiles[50]) / 2;
    }

    public static double getTpmPercentile(final Map<String,double[]> transPercentilesMap, final String transName, double tpm)
    {
        final double[] transPercentiles = transPercentilesMap.get(transName);

        if(transPercentiles == null)
            return -1;

        if(tpm < transPercentiles[0])
            return 0;
        else if(tpm > transPercentiles[transPercentiles.length - 1])
            return transPercentiles.length - 1;

        for(int i = 0; i < transPercentiles.length - 1; ++i)
        {
            if(tpm >= transPercentiles[i] && tpm <= transPercentiles[i + 1])
            {
                if(transPercentiles[i + 1] == transPercentiles[i])
                    return i;

                double upperFactor = (tpm - transPercentiles[i]) / (transPercentiles[i + 1] - transPercentiles[i]);
                return upperFactor * (i + 1) + (1 - upperFactor) * i;
            }
        }

        return -1;
    }

}
