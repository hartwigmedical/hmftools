package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.SAMPLE_GENE_PERCENTILES;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.AnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class SampleGenePercentiles
{
    private final CohortConfig mConfig;

    private final Map<String,Map<String,double[]>> mCancerTypesGeneDistribution;

    private BufferedWriter mWriter;

    private static final String PAN_CANCER = "ALL";
    private static final String UNKNOWN_CANCER = "Unknown";
    private static final int FLD_CANCER_TYPE = 0;
    private static final int FLD_FILENAME = 1;

    public SampleGenePercentiles(final CohortConfig config)
    {
        mConfig = config;
        mCancerTypesGeneDistribution = Maps.newHashMap();

        loadCancerPercentileData();

        mWriter = null;
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

                double[] percentileData = new double[PERCENTILE_COUNT];

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

    private void loadCancerPercentileData()
    {
        if(mConfig.Expression.CancerGeneFiles == null)
            return;

        final Map<String, String> cancerGeneFilenames = loadCancerGeneDistributionFilenames(mConfig.Expression.CancerGeneFiles);

        for(Map.Entry<String, String> entry : cancerGeneFilenames.entrySet())
        {
            final String cancerType = entry.getKey();
            final String filename = entry.getValue();

            ISF_LOGGER.debug("cancerType({}) loading gene records", cancerType);

            final Map<String, double[]> geneDistributionMap = Maps.newHashMap();
            loadCohortDistribution(filename, geneDistributionMap, "gene", PERCENTILE_COUNT + 2, mConfig.RestrictedGeneIds);

            ISF_LOGGER.info("cancerType({}) loaded {} gene records", cancerType, geneDistributionMap.size());

            mCancerTypesGeneDistribution.put(cancerType, geneDistributionMap);
        }
    }

    public void processSampleFiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, SAMPLE_GENE_PERCENTILES, filenames))
            return;

        initialiseWriter();

        ISF_LOGGER.info("processing {} samples gene files", mConfig.SampleData.SampleIds.size());

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path genesFile = filenames.get(i);

            processSampleFile(sampleId, genesFile);
            ISF_LOGGER.debug("{}: sample({}) processed gene file", i, sampleId);
        }

        ISF_LOGGER.info("processed {} samples gene files", mConfig.SampleData.SampleIds.size());

        closeBufferedWriter(mWriter);
    }

    private void processSampleFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            final String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);

            final Map<String,double[]> panCancerPercentilesMap = getCancerTypePercentilesMap(PAN_CANCER);
            final Map<String,double[]> cancerPercentilesMap = getCancerTypePercentilesMap(cancerType);

            final List<String> sampleGeneIds = mConfig.SampleData.SampleGeneIds.get(sampleId);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMap.get(FLD_GENE_NAME);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                if(sampleGeneIds != null && !sampleGeneIds.contains(geneId))
                    continue;

                double tpm = Double.parseDouble(items[tpmIndex]);

                if(tpm < mConfig.Expression.TpmThreshold)
                    continue;

                final String geneName = items[geneNameIndex];

                double panCancerPerc = 0;
                double cancerPerc = 0;
                double panCancerMedian = 0;
                double cancerMedian = 0;

                if(panCancerPercentilesMap != null)
                {
                    panCancerPerc = getTpmPercentile(panCancerPercentilesMap, geneId, tpm);
                    panCancerMedian = getTpmMedian(panCancerPercentilesMap, geneId);
                }

                if(cancerPercentilesMap != null)
                {
                    cancerPerc = getTpmPercentile(cancerPercentilesMap, geneId, tpm);
                    cancerMedian = getTpmMedian(cancerPercentilesMap, geneId);
                }

                writeSamplePercentileData(
                        sampleId, cancerType, geneId, geneName, tpm, panCancerPerc, cancerPerc, panCancerMedian, cancerMedian);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
            return;
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

        return getPercentile(transPercentiles, tpm);
    }

    private Map<String,double[]> getCancerTypePercentilesMap(final String cancerType)
    {
        final Map<String,double[]> transPercentilesMap = mCancerTypesGeneDistribution.get(cancerType);

        if(transPercentilesMap != null)
            return transPercentilesMap;

        if(!mCancerTypesGeneDistribution.containsKey(UNKNOWN_CANCER))
            return null;

        return mCancerTypesGeneDistribution.get(UNKNOWN_CANCER);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("sample_gene_perc_data.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,CancerType,GeneId,GeneName,TPM");

            if(!mCancerTypesGeneDistribution.isEmpty())
                mWriter.write(",CohortMedianTPM,CancerMedianTPM,CohortPercentile,CancerPercentile");

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void writeSamplePercentileData(
            final String sampleId, final String cancerType, final String geneId, final String geneName, double tpm,
            double panCancerPerc, double cancerPerc, double panCancerMedian, double cancerMedian)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%6.3e", sampleId, cancerType, geneId, geneName, tpm));

            if(!mCancerTypesGeneDistribution.isEmpty())
                mWriter.write(String.format(",%6.3e,%6.3e,%.1f,%.1f", panCancerMedian, cancerMedian, panCancerPerc, cancerPerc));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write sample gene percentile data file: {}", e.toString());
        }
    }

    private Map<String,String> loadCancerGeneDistributionFilenames(final String inputFile)
    {
        final Map<String,String> cancerGeneFileMap = Maps.newHashMap();

        if(!Files.exists(Paths.get(inputFile)))
        {
            ISF_LOGGER.error("invalid cancer-gene distribution file list");
            return cancerGeneFileMap;
        }

        try
        {
            final List<String> lines = Files.readAllLines(new File(inputFile).toPath());

            for(String data : lines)
            {
                final String[] items = data.split(",");
                final String cancerType = items[FLD_CANCER_TYPE];
                final String filename = items[FLD_FILENAME];

                cancerGeneFileMap.put(cancerType, filename);
            }

            ISF_LOGGER.info("loaded {} cancer-gene distribution filenames from file({})", cancerGeneFileMap.size(), inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load cancer-gene distribution file({}): {}", inputFile, e.toString());
        }

        return cancerGeneFileMap;
    }

}
