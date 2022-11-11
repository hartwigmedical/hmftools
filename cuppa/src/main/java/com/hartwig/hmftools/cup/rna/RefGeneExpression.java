package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_CANCER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_SAMPLE;
import static com.hartwig.hmftools.common.cuppa.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.GENE_EXP_IGNORE_FIELDS;
import static com.hartwig.hmftools.cup.rna.GeneExpressionDataLoader.loadGeneExpressionMatrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefClassifier;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class RefGeneExpression implements RefClassifier
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private Matrix mCancerGeneExpression;
    private Matrix mSampleGeneExpression;
    private final Map<String,Integer> mSampleTpmIndex;
    private final List<String> mSampleNames;
    private final List<String> mGeneIds;
    private final List<String> mGeneNames;
    private final List<String> mCancerTypes;

    private final boolean mTpmInLogForm;

    private static final String TPM_IN_LOG_FORM = "tpm_as_log";

    public RefGeneExpression(final RefDataConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mCancerGeneExpression = null;
        mSampleGeneExpression = null;
        mSampleTpmIndex = Maps.newHashMap();
        mGeneIds = Lists.newArrayList();
        mGeneNames = Lists.newArrayList();
        mCancerTypes = Lists.newArrayList();
        mSampleNames = Lists.newArrayList();

        mTpmInLogForm = cmd.hasOption(TPM_IN_LOG_FORM);
    }

    public static void addCmdLineArgs(@NotNull Options options)
    {
        options.addOption(TPM_IN_LOG_FORM, false, "Expect TPM in log form");
    }

    public CategoryType categoryType() { return GENE_EXP; }

    public static boolean requiresBuild(final RefDataConfig config)
    {
        return config.Categories.contains(GENE_EXP) || !config.GeneExpMatrixFile.isEmpty();
    }

    public void buildRefDataSets()
    {
        CUP_LOGGER.debug("loading RNA gene expression data");

        // loadRefRnaGeneExpression(mConfig.GeneExpMatrixFile);

        Matrix sampleGeneExpression = loadGeneExpressionMatrix(
                mConfig.GeneExpMatrixFile, mSampleTpmIndex, mSampleNames, mGeneIds, mGeneNames);

        // CUP_LOGGER.debug("loaded {} gene for expression ref data", geneIds.size());

        // limit to the ref samples and then creating columns in the same order as loaded
        mSampleGeneExpression = new Matrix(mSampleDataCache.RefSampleDataList.size(), mGeneIds.size());

        if(mSampleGeneExpression == null)
        {
            CUP_LOGGER.warn("RNA gene expression data load failed");
            return;
        }

        double[][] sampleData = mSampleGeneExpression.getData();

        for(int i = 0; i < mSampleDataCache.RefSampleDataList.size(); ++i)
        {
            SampleData refSample = mSampleDataCache.RefSampleDataList.get(i);

            if(!mSampleNames.contains(refSample.Id))
            {
                CUP_LOGGER.error("sample({}) missing from gene expression matrix", refSample.Id);
                continue;
            }

            int countsIndex = mSampleTpmIndex.get(refSample.Id);
            double[] sampleTPMs = sampleGeneExpression.getRow(countsIndex);

            for(int b = 0; b < sampleTPMs.length; ++b)
            {
                // transformation: double logTpm = log(adjTpm + 1);
                double value = sampleTPMs[b];

                if(mTpmInLogForm)
                    sampleData[i][b] = exp(value) - 1;
                else
                    sampleData[i][b] = value;
            }
        }

        if(mConfig.NoiseAdjustments.hasNoiseAllocation(EXPRESSION_PAIRWISE))
        {
            int noiseAllocation = mConfig.NoiseAdjustments.getNoiseAllocation(EXPRESSION_PAIRWISE);

            final double[] geneExpMedians = NoiseRefCache.generateMedianValues(mSampleGeneExpression);
            double medianTotals = sumVector(geneExpMedians);

            CUP_LOGGER.debug("applying noise({}) to gene expression sample TPMs, medianTotal({})",
                    noiseAllocation, format("%.0f", medianTotals));

            NoiseRefCache.applyNoise(mSampleGeneExpression, geneExpMedians, noiseAllocation);

            // write medians so they can be applied to a new sample if required
            mConfig.NoiseAdjustments.addNoiseData(EXPRESSION_PAIRWISE, geneExpMedians);
        }

        if(mTpmInLogForm)
        {
            // convert back
            final double[][] data = mSampleGeneExpression.getData();
            for(int s = 0; s < mSampleGeneExpression.Rows; ++s)
            {
                for(int b = 0; b < mSampleGeneExpression.Cols; ++b)
                {
                    double tpm = data[s][b];
                    data[s][b] = log(tpm + 1);
                }
            }
        }

        CUP_LOGGER.debug("writing RNA gene expression sample data");
        writeMatrix(mSampleGeneExpression, mSampleNames, REF_FILE_GENE_EXP_SAMPLE);

        // create a per-cancer type sum of sample counts
        mCancerGeneExpression = new Matrix(mSampleDataCache.RefCancerSampleData.size(), mSampleGeneExpression.Cols);
        final double[][] cancerMatrixData = mCancerGeneExpression.getData();

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

                final double[] sampleTPMs = mSampleGeneExpression.getRow(countsIndex);

                for(int b = 0; b < sampleTPMs.length; ++b)
                {
                    cancerMatrixData[cancerIndex][b] += sampleTPMs[b];
                }
            }
        }

        CUP_LOGGER.debug("writing RNA gene expression cohort data");

        writeMatrix(mCancerGeneExpression, mCancerTypes, REF_FILE_GENE_EXP_CANCER);
    }

    private void loadRefRnaGeneExpression(final String filename)
    {
        if(filename.isEmpty())
            return;

        Matrix sampleGeneExpression = loadMatrixDataFile(filename, mSampleTpmIndex, GENE_EXP_IGNORE_FIELDS, true);

        // first populate the gene info and sample names
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();

            String[] columns = header.split(DATA_DELIM, -1);

            if(columns.length < 3 || !columns[0].equals(FLD_GENE_ID) || !columns[1].equals(FLD_GENE_NAME))
            {
                CUP_LOGGER.error("invalid gene expression file header");
                return;
            }

            // assumes GeneId, GeneName, Sample1, Sample2 etc
            for(int i = 2; i < columns.length; ++i)
            {
                mSampleNames.add(columns[i]);
            }

            String line = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DATA_DELIM, -1);
                mGeneIds.add(items[0]);
                mGeneNames.add(items[1]);
            }

            CUP_LOGGER.debug("loaded {} gene for expression ref data", mGeneIds.size());

            // make a matrix by first checking the ref samples and then creating columns in the same order as loaded
            mSampleGeneExpression = new Matrix(mSampleDataCache.RefSampleDataList.size(), mGeneIds.size());
            double[][] sampleData = mSampleGeneExpression.getData();

            for(int i = 0; i < mSampleDataCache.RefSampleDataList.size(); ++i)
            {
                SampleData refSample = mSampleDataCache.RefSampleDataList.get(i);

                if(!mSampleNames.contains(refSample.Id))
                {
                    CUP_LOGGER.error("sample({}) missing from gene expression matrix", refSample.Id);
                    continue;
                }

                int countsIndex = mSampleTpmIndex.get(refSample.Id);
                double[] sampleTPMs = sampleGeneExpression.getRow(countsIndex);

                for(int b = 0; b < sampleTPMs.length; ++b)
                {
                    // transformation: double logTpm = log(adjTpm + 1);
                    double value = sampleTPMs[b];

                    if(mTpmInLogForm)
                        sampleData[i][b] = exp(value) - 1;
                    else
                        sampleData[i][b] = value;
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.debug("failed to load RNA expression ref data from {}: {}", filename, e.toString());
        }
    }

    private void writeMatrix(final Matrix tpmMatrix, final List<String> headers, final String fileId)
    {
        try
        {
            final String filename = mConfig.OutputDir + fileId;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName");

            for(final String header : headers)
            {
                writer.write(format(",%s", header));
            }

            writer.newLine();

            // samples / cancer types are in rows, TPMs in columns but write transposed
            final double[][] matrixData = tpmMatrix.getData();

            for(int b = 0; b < tpmMatrix.Cols; ++b)
            {
                writer.write(format("%s,%s", mGeneIds.get(b), mGeneNames.get(b)));

                for(int i = 0; i < tpmMatrix.Rows; ++i)
                {
                    double tpm = matrixData[i][b];

                    // write decimal is most efficient form
                    if(tpm == 0)
                        writer.write(",");
                    else if(tpm > 999 || tpm < 0.001)
                        writer.write(format(",%6.3e", tpm));
                    else
                        writer.write(format(",%.4g", tpm));
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
}
