package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class ExpressionCohortDistribution
{
    // for re-writing expression matrix as pan-cancer and per-cancer median values and percentiles
    private final CohortConfig mConfig;
    private BufferedWriter mWriter;

    private Matrix mExpressionMatrix;
    private final Map<String,Integer> mSampleIndexMap; // index of sample into matrix columns
    private final Map<String,Integer> mGeneTransIdIndexMap;  // index of geneId into matrix rows
    private final List<String> mGeneTransIds;
    private final Map<String,String> mIdNameMap;

    private boolean mTranscriptScope;
    private boolean mWritePerCancer;
    private boolean mWritePercentiles;

    private static final String DIST_BY_CANCER_TYPE = "exp_by_cancer";
    private static final String TRANSCRIPT_SCOPE = "exp_by_transcript";
    private static final String WRITE_PERCENTILES = "exp_percentiles";

    public ExpressionCohortDistribution(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mExpressionMatrix = null;
        mSampleIndexMap = Maps.newHashMap();
        mGeneTransIdIndexMap = Maps.newHashMap();
        mGeneTransIds = Lists.newArrayList();
        mIdNameMap = Maps.newHashMap();

        mTranscriptScope = cmd.hasOption(TRANSCRIPT_SCOPE);
        mWritePercentiles = cmd.hasOption(WRITE_PERCENTILES);
        mWritePerCancer = !mWritePercentiles || cmd.hasOption(DIST_BY_CANCER_TYPE);

        mWriter = null;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(DIST_BY_CANCER_TYPE, false, "Produce per-cancer expression distributions");
        options.addOption(TRANSCRIPT_SCOPE, false, "Produce per-cancer expression distributions");
        options.addOption(WRITE_PERCENTILES, false, "Write expression percentiles");
    }

    public void produceCohortData()
    {
        if(!loadExpressionMatrix())
            return;

        initialiseWriter();

        int nextLog = 1000;

        for(int i = 0; i < mGeneTransIds.size(); ++i)
        {
            produceCohortData(i);

            if(i >= nextLog)
            {
                nextLog += 1000;
                ISF_LOGGER.info("processed {} {}", i, mTranscriptScope ? "transcripts" : "genes");
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void produceCohortData(int index)
    {
        final String geneTransId = mGeneTransIds.get(index);

        final String geneTransName = mIdNameMap.get(geneTransId);
        int matrixIndex = mGeneTransIdIndexMap.get(geneTransId);

        final double[] tpmValues = mExpressionMatrix.getRow(matrixIndex);

        try
        {
            // first write out pan-cancer medians and optionally percentiles
            final List<Double> allValues = Lists.newArrayListWithExpectedSize(mExpressionMatrix.Cols);

            for(int s = 0; s < tpmValues.length; ++s)
            {
                addSortedTpm(allValues, tpmValues[s]);
            }

            if(mWritePercentiles)
            {
                writeCancerValues(geneTransId, geneTransName, "ALL", allValues);
            }
            else
            {
                double allMedian = calculatedMedian(allValues);

                if(mTranscriptScope)
                    mWriter.write(String.format("%s,%s", geneTransName, geneTransId));
                else
                    mWriter.write(String.format("%s,%s", geneTransId, geneTransName));

                mWriter.write(String.format(",%6.3e", allMedian));
            }

            if(mWritePerCancer)
            {
                for(Map.Entry<String, List<String>> entry : mConfig.SampleData.CancerTypeSamples.entrySet())
                {
                    final String cancerType = entry.getKey();
                    final List<String> samples = entry.getValue();

                    final List<Double> cancerValues = Lists.newArrayListWithExpectedSize(samples.size());

                    for(final String sampleId : samples)
                    {
                        Integer sampleIndex = mSampleIndexMap.get(sampleId);

                        if(sampleIndex == null)
                            continue;

                        addSortedTpm(cancerValues, tpmValues[sampleIndex]);
                    }

                    if(mWritePercentiles)
                    {
                        writeCancerValues(geneTransId, geneTransName, cancerType, cancerValues);
                    }
                    else
                    {
                        double cancerMedian = calculatedMedian(cancerValues);
                        mWriter.write(String.format(",%6.3e", cancerMedian));
                    }
                }
            }

            if(!mWritePercentiles)
                mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private double calculatedMedian(final List<Double> values)
    {
        if((values.size() % 2) == 0)
        {
            int medianIndexLow = values.size() / 2;
            int medianIndexHigh = medianIndexLow + 1;
            return (values.get(medianIndexLow) + values.get(medianIndexHigh)) * 0.5;
        }
        else
        {
            int medianIndex = values.size() / 2;
            return values.get(medianIndex);
        }
    }

    private void addSortedTpm(final List<Double> values, double tpm)
    {
        int index = 0;
        while(index < values.size())
        {
            if(tpm < values.get(index))
                break;

            ++index;
        }

        values.add(index, tpm);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mTranscriptScope ?
                    mConfig.formCohortFilename("transcript_distribution.csv") : mConfig.formCohortFilename("gene_distribution.csv");

            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("GeneId,GeneName");

            if(mTranscriptScope)
                mWriter.write(",TransName");

            if(mWritePercentiles)
            {
                if(mWritePerCancer)
                    mWriter.write(",CancerType");

                mWriter.write(",Median");

                for(int i = 0; i < PERCENTILE_COUNT; ++i)
                {
                    mWriter.write(String.format(",Pct_%d", i));
                }
            }
            else
            {
                //  cancer types for the remaining columns
                mWriter.write(",All");

                for(String cancerType : mConfig.SampleData.CancerTypeSamples.keySet())
                {
                    mWriter.write(String.format(",%s", cancerType));
                }
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private void writeCancerValues(final String geneTransId, final String geneTransName, final String cancerType, final List<Double> tpmList)
    {
        try
        {
            if(mTranscriptScope)
                mWriter.write(String.format("%s,%s", geneTransName, geneTransId));
            else
                mWriter.write(String.format("%s,%s", geneTransId, geneTransName));

            if(mWritePerCancer)
                mWriter.write(String.format(",%s", cancerType));

            final double[] tmpArray = convertList(tpmList);

            final double[] percentileValues = new double[PERCENTILE_COUNT];
            calcPercentileValues(tmpArray, percentileValues);

            double medianValue = calculatedMedian(tpmList);

            mWriter.write(String.format(",%6.3e", medianValue));

            for (int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                mWriter.write(String.format(",%6.3e", percentileValues[i]));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private boolean loadExpressionMatrix()
    {
        // keep track of gene/transcript ids and names
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.Expression.GeneExpMatrixFile));

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMapIndex = createFieldsIndexMap(header, DELIMITER);

            int geneIdIndex = fieldsMapIndex.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMapIndex.get(FLD_GENE_NAME);
            int transNameIndex = fieldsMapIndex.containsKey(FLD_TRANS_NAME) ? fieldsMapIndex.get(FLD_TRANS_NAME) : -1;
            mTranscriptScope = transNameIndex >= 0;

            String line = fileReader.readLine();
            int geneRowIndex = 0;

            while(line != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneId = items[geneIdIndex];

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                {
                    geneRowIndex++;
                    line = fileReader.readLine();
                    continue;
                }

                final String geneName = items[geneNameIndex];

                if(mTranscriptScope)
                {
                    final String transName = items[transNameIndex];
                    mGeneTransIds.add(transName);
                    mGeneTransIdIndexMap.put(transName, geneRowIndex);
                    mIdNameMap.put(transName, String.format("%s,%s", geneId, geneName));
                }
                else
                {
                    mGeneTransIds.add(geneId);
                    mGeneTransIdIndexMap.put(geneId, geneRowIndex);
                    mIdNameMap.put(geneId, geneName);
                }

                geneRowIndex++;
                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.debug("failed to load RNA expression data from {}: {}", mConfig.Expression.GeneExpMatrixFile, e.toString());
            return false;
        }

        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME);

        if(mTranscriptScope)
            ignoreFields.add(FLD_TRANS_NAME);

        mExpressionMatrix = loadMatrixDataFile(mConfig.Expression.GeneExpMatrixFile, mSampleIndexMap, ignoreFields);

        ISF_LOGGER.debug("loaded genes({}) and {} samples expression matrix data", mGeneTransIds.size(), mExpressionMatrix.Cols);
        return true;
    }

}
