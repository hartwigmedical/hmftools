package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class RnaExpressionMatrix
{
    private final Matrix mExpressionMatrix;
    private final Map<String,Integer> mSampleIdIndexMap; // sampleId to column index
    private final Map<String,Integer> mGeneTransIdIndexMap; // geneId or transName to matrix row index
    private final boolean mTranscriptScope;
    private double mMinTpmThreshold;

    public static final String EXPRESSION_SCOPE_GENE = "GENE";
    public static final String EXPRESSION_SCOPE_TRANS = "TRANS";

    public RnaExpressionMatrix(final String expressionFile, final String scope)
    {
        mTranscriptScope = scope.equals(EXPRESSION_SCOPE_TRANS);
        mMinTpmThreshold = 0;

        mSampleIdIndexMap = Maps.newHashMap();
        mGeneTransIdIndexMap = Maps.newHashMap();
        mExpressionMatrix = load(expressionFile);
    }

    public void setMinTpmThreshold(double threshold) { mMinTpmThreshold = threshold; }

    public boolean hasValidData()
    {
        return mExpressionMatrix != null
                && mGeneTransIdIndexMap.size() == mExpressionMatrix.Rows
                && mSampleIdIndexMap.size() == mExpressionMatrix.Cols;
    }

    public static final double INVALID_EXP = -1;

    public double getExpression(final String geneTransId, final String sampleId)
    {
        Integer geneIndex = mGeneTransIdIndexMap.get(geneTransId);
        Integer sampleIndex = mSampleIdIndexMap.get(sampleId);

        if(geneIndex == null || sampleIndex == null)
            return INVALID_EXP;

        return mExpressionMatrix.get(geneIndex, sampleIndex);
    }

    public boolean hasGeneId(final String geneId) { return mGeneTransIdIndexMap.containsKey(geneId); }
    public boolean hasTransName(final String transName) { return mGeneTransIdIndexMap.containsKey(transName); }
    public boolean hasSampleId(final String sampleId) { return mSampleIdIndexMap.containsKey(sampleId); }

    private String scopeType() { return mTranscriptScope ? "transcript" : "gene"; }

    private Matrix load(final String filename)
    {
        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME);

        if(mTranscriptScope)
            ignoreFields.add(FLD_TRANS_NAME);

        Matrix matrix = loadMatrixDataFile(filename, mSampleIdIndexMap, ignoreFields);

        List<Integer> excludedRows = Lists.newArrayList();

        if(mMinTpmThreshold > 0)
            matrix = checkExcludedRows(matrix, excludedRows);

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileData.remove(0);

            int idCol = mTranscriptScope ? fieldsIndexMap.get(FLD_TRANS_NAME) : fieldsIndexMap.get(FLD_GENE_ID);

            for(int i = 0; i < fileData.size(); ++i)
            {
                final String[] items = fileData.get(i).split(DELIMITER, -1);
                final String idValue = items[idCol];
                mGeneTransIdIndexMap.put(idValue, i);
            }
        }
        catch (IOException e)
        {
            RNA_LOGGER.error("failed to read RNA {} expression data from file({}): {}", scopeType(), filename, e.toString());
            return null;
        }

        RNA_LOGGER.info("loaded expression for {} {}s and {} samples from file({})",
                mGeneTransIdIndexMap.size(), scopeType(), mSampleIdIndexMap.size(), filename);

        return matrix;
    }

    private Matrix checkExcludedRows(final Matrix matrix, final List<Integer> excludedRows)
    {
        // filter out transcrpts where no TPM exceeds the threshold
        final double[][] data = matrix.getData();

        for(int row = 0; row < matrix.Rows; ++row)
        {
            boolean allLow = true;

            for(int col = 0; col < matrix.Cols; ++col)
            {
                if(data[row][col] > mMinTpmThreshold)
                {
                    allLow = false;
                    break;
                }
            }

            if(allLow)
                excludedRows.add(row);
        }

        if(excludedRows.isEmpty())
            return matrix;

        RNA_LOGGER.debug("excluding {} {}(s) with low TPM vs threshold({})", excludedRows.size(), scopeType(), mMinTpmThreshold);

        int newRowCount = matrix.Rows - excludedRows.size();
        Matrix newMatrix = new Matrix(newRowCount, matrix.Cols);
        final double[][] newData = matrix.getData();

        for(int row = 0; row < matrix.Rows; ++row)
        {
            if(excludedRows.contains(row))
                continue;

            for(int col = 0; col < matrix.Rows; ++col)
            {
                newData[row][col] = data[row][col];
            }
        }

        return newMatrix;
    }
}
