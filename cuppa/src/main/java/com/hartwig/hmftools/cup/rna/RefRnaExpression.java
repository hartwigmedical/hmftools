package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_CANCER;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class RefRnaExpression implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private Matrix mGeneCancerExpressionData;
    private Matrix mGeneSampleExpressionData;
    private final Map<String,Integer> mSampleTpmIndex;
    private final List<String> mGeneIds;
    private final List<String> mGeneNames;
    private final List<String> mSampleIds;
    private final List<String> mCancerTypes;
    private final double mTpmLogCutoff;

    public static final String TPM_LOG_CUTOFF = "rna_tpm_log_cutoff";

    public RefRnaExpression(final RefDataConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mGeneCancerExpressionData = null;
        mGeneSampleExpressionData = null;
        mSampleTpmIndex = Maps.newHashMap();
        mGeneIds = Lists.newArrayList();
        mGeneNames = Lists.newArrayList();
        mSampleIds = Lists.newArrayList();
        mCancerTypes = Lists.newArrayList();

        mTpmLogCutoff = Double.parseDouble(cmd.getOptionValue(TPM_LOG_CUTOFF, "0"));
    }

    public static void addCmdLineArgs(@NotNull Options options)
    {
        options.addOption(TPM_LOG_CUTOFF, true, "RNA TPM cut-off in log scale (default=0, not applied)");
    }

    public CategoryType categoryType() { return GENE_EXP; }

    public static boolean requiresBuild(final RefDataConfig config) { return !config.RefRnaGeneExpFile.isEmpty(); }

    public void buildRefDataSets()
    {
        CUP_LOGGER.debug("loading RNA gene expression data");

        loadRefRnaGeneExpression(mConfig.RefRnaGeneExpFile);

        if(mGeneSampleExpressionData == null)
        {
            CUP_LOGGER.warn("RNA gene expression data load failed");
            return;
        }

        mGeneCancerExpressionData = new Matrix(mGeneSampleExpressionData.Rows, mSampleDataCache.RefCancerSampleData.size());
        final double[][] cancerMatrixData = mGeneCancerExpressionData.getData();

        for(Map.Entry<String,List<SampleData>> entry : mSampleDataCache.RefCancerSampleData.entrySet())
        {
            final String refCancerType = entry.getKey();
            mCancerTypes.add(refCancerType);
            int cancerIndex = mCancerTypes.size() - 1;

            for(final SampleData sample : entry.getValue())
            {
                Integer countsIndex = mSampleTpmIndex.get(sample.Id);

                if(countsIndex == null)
                {
                    CUP_LOGGER.warn("sample({}) missing gene expression data", sample.Id);
                    continue;
                }

                final double[] sampleTPMs = mGeneSampleExpressionData.getCol(countsIndex);

                for(int b = 0; b < sampleTPMs.length; ++b)
                {
                    cancerMatrixData[b][cancerIndex] += sampleTPMs[b];
                }
            }
        }

        CUP_LOGGER.debug("writing RNA gene expression reference data");

        writeMatrixData(mGeneCancerExpressionData, mCancerTypes, REF_FILE_GENE_EXP_CANCER);
    }

    private void loadRefRnaGeneExpression(final String filename)
    {
        if(filename.isEmpty())
            return;

        final List<String> ignoreFields = Lists.newArrayList("GeneId", "GeneName");

        mGeneSampleExpressionData = loadMatrixDataFile(filename, mSampleTpmIndex, ignoreFields);
        mGeneSampleExpressionData.cacheTranspose();

        // populate the gene info
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMapIndex = createFieldsIndexMap(header, DATA_DELIM);

            int geneIdIndex = fieldsMapIndex.get("GeneId");
            int geneNameIndex = fieldsMapIndex.get("GeneName");

            String line = fileReader.readLine();

            while(line != null)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                mGeneIds.add(items[geneIdIndex]);
                mGeneNames.add(items[geneNameIndex]);

                line = fileReader.readLine();
            }

            CUP_LOGGER.debug("loaded {} gene for expression ref data", mGeneIds.size());
        }
        catch (IOException e)
        {
            CUP_LOGGER.debug("failed to load RNA expression ref data from {}: {}", filename, e.toString());
        }
    }

    private void writeMatrixData(final Matrix tpmMatrix, final List<String> headers, final String fileId)
    {
        try
        {
            final String filename = mConfig.OutputDir + fileId;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName");

            for(final String header : headers)
            {
                writer.write(String.format(",%s", header));
            }

            writer.newLine();

            final double[][] matrixData = tpmMatrix.getData();

            for(int i = 0; i < tpmMatrix.Rows; ++i)
            {
                writer.write(String.format("%s,%s", mGeneIds.get(i), mGeneNames.get(i)));

                for(int j = 0; j < tpmMatrix.Cols; ++j)
                {
                    writer.write(String.format(",%.2f", matrixData[i][j]));
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref RNA gene expression output: {}", e.toString());
        }
    }

    public static boolean loadRefPercentileData(
            final String filename, final Map<String,Map<String,double[]>> refPercentiles, final Map<String,String> geneIdNameMap)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();

            // SampleId,CancerType,GeneId,GeneName,TPM
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            int cancerCol = fieldsIndexMap.get("CancerType");
            int geneIdCol = fieldsIndexMap.get("GeneId");
            int geneNameCol = fieldsIndexMap.get("GeneName");

            String line = fileReader.readLine();

            String currentGeneId = "";
            Map<String,double[]> cancerPercMap = null;

            while(line != null)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String geneId = items[geneIdCol];
                final String geneName = items[geneNameCol];

                if(!currentGeneId.equals(geneId))
                {
                    currentGeneId = geneId;
                    geneIdNameMap.put(geneId,geneName);
                    cancerPercMap = Maps.newHashMap();
                    refPercentiles.put(geneId,cancerPercMap);
                }

                final String cancerType = items[cancerCol];
                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 3;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                cancerPercMap.put(cancerType, percentileData);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA exp cancer percentile data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }
}
