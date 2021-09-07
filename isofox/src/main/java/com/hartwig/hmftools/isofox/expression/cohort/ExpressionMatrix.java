package com.hartwig.hmftools.isofox.expression.cohort;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.TRANSCRIPT_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.isofox.cohort.AnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class ExpressionMatrix
{
    private final AnalysisType mType;
    private final CohortConfig mConfig;
    
    private Matrix mExpressionMatrix;
    private final List<String> mGeneIds;
    private final List<String> mGeneNames;
    private final List<String> mTranscriptNames;

    public ExpressionMatrix(final CohortConfig config, final AnalysisType type)
    {
        mType = type;
        mConfig = config;

        mExpressionMatrix = null;
        mGeneIds = Lists.newArrayList();
        mGeneNames = Lists.newArrayList();
        mTranscriptNames = Lists.newArrayList();
    }

    public void processSamples()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, mType, filenames))
            return;

        final String typeStr = mType == GENE_EXPRESSION_MATRIX ? "gene" : "transcript";

        ISF_LOGGER.info("processing {} samples {} files", mConfig.SampleData.SampleIds.size(), typeStr);

        Map<String,Integer> fieldsMap = Maps.newHashMap();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path sampleFile = filenames.get(i);

            processSampleFile(i, sampleFile, fieldsMap);
            ISF_LOGGER.debug("{}: sample({}) processed {} file", i, sampleId, typeStr);
        }

        ISF_LOGGER.info("processed {} samples {} files", mConfig.SampleData.SampleIds.size(), typeStr);

        writeMatrixData();
    }

    private void processSampleFile(int sampleIndex, final Path filename, final Map<String,Integer> fieldsMap)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(fieldsMap.isEmpty())
                fieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            boolean isTranscriptLevel = mType == TRANSCRIPT_EXPRESSION_MATRIX;

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMap.get(FLD_GENE_NAME);

            int transNameIndex = isTranscriptLevel ? fieldsMap.get(FLD_TRANS_NAME) : -1;
            int tpmIndex = fieldsMap.get(FLD_TPM);

            int itemIndex = 0;
            List<String> itemCache = isTranscriptLevel ? mTranscriptNames : mGeneIds;
            int manualLookupCount = 0;

            double[][] sampleMatrixData = null;

            if(mExpressionMatrix == null)
            {
                int samplesCount = mConfig.SampleData.SampleIds.size();
                int expressionItemCount = lines.size();

                // cull rows based on any restrictions in place
                if(!mConfig.RestrictedGeneIds.isEmpty())
                {
                    long itemCount = lines.stream().map(x -> x.split(DELIMITER, -1)[geneIdIndex])
                            .filter(x -> mConfig.RestrictedGeneIds.contains(x))
                            .count();

                    expressionItemCount = (int)itemCount;
                }

                ISF_LOGGER.debug("building gene expression matrix: genes({}) samples({})", expressionItemCount, samplesCount);

                mExpressionMatrix = new Matrix(expressionItemCount, samplesCount);
            }

            sampleMatrixData = mExpressionMatrix.getData();

            boolean buildIndexes = itemCache.isEmpty();

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];
                final String transName = transNameIndex >= 0 ? items[transNameIndex] : "";

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                if(buildIndexes)
                {
                    mGeneIds.add(itemIndex, geneId);
                    mGeneNames.add(itemIndex, items[geneNameIndex]);

                    if(transNameIndex >= 0)
                    {
                        mTranscriptNames.add(itemIndex, transName);
                    }
                }
                else
                {
                    final String itemId = isTranscriptLevel ? transName : geneId;

                    if(itemIndex >= itemCache.size() || !itemCache.get(itemIndex).equals(itemId))
                    {
                        ++manualLookupCount;

                        // locate manually
                        boolean found = false;
                        for(itemIndex = 0; itemIndex < itemCache.size(); ++itemIndex)
                        {
                            if(itemCache.get(itemIndex).equals(itemId))
                            {
                                found = true;
                                break;
                            }
                        }

                        if(!found)
                        {
                            ISF_LOGGER.error("item({}) not present in item cache samples", itemId);
                            return;
                        }
                    }
                }

                double tpm = Double.parseDouble(items[tpmIndex]);

                if(mConfig.Expression.TpmThreshold > 0 && tpm < mConfig.Expression.TpmThreshold)
                {
                    ++itemIndex;
                    continue;
                }

                if(mConfig.Expression.UseLogTpm)
                    tpm = log(tpm + 1);

                sampleMatrixData[itemIndex][sampleIndex] = tpm;
                ++itemIndex;
            }

            if(manualLookupCount > 0)
            {
                ISF_LOGGER.debug("required {} manual gene look-ups", manualLookupCount);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private void writeMatrixData()
    {
        try
        {
            final String fileType = (mType == TRANSCRIPT_EXPRESSION_MATRIX ? "transcript" : "gene") + "_expression_matrix.csv";
            final String filename = mConfig.formCohortFilename(fileType);

            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName");

            if(mType == TRANSCRIPT_EXPRESSION_MATRIX)
                writer.write(",TransName");

            for(final String sampleId : mConfig.SampleData.SampleIds)
            {
                writer.write(String.format(",%s", sampleId));
            }

            writer.newLine();

            final double[][] matrixData = mExpressionMatrix.getData();

            for(int row = 0; row < mExpressionMatrix.Rows; ++row)
            {
                if(mConfig.Expression.ApplyTpmWriteLimit)
                {
                    // skip a row if all sample entries are below the specified threshold
                    boolean allLow = true;

                    for(int col = 0; col < mExpressionMatrix.Cols; ++col)
                    {
                        if(matrixData[row][col] >= mConfig.Expression.TpmThreshold)
                        {
                            allLow = false;
                            break;
                        }
                    }

                    if(allLow)
                        continue;
                }

                writer.write(String.format("%s,%s", mGeneIds.get(row), mGeneNames.get(row)));

                if(mType == TRANSCRIPT_EXPRESSION_MATRIX)
                    writer.write(String.format(",%s", mTranscriptNames.get(row)));

                for(int j = 0; j < mExpressionMatrix.Cols; ++j)
                {
                    // write decimal is most efficient form
                    if(matrixData[row][j] == 0)
                        writer.write(",0");
                    else if(matrixData[row][j] > 999 || matrixData[row][j] < 0.001)
                        writer.write(String.format(",%6.3e", matrixData[row][j]));
                    else
                        writer.write(String.format(",%.4g", matrixData[row][j]));
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write expression matrix output: {}", e.toString());
        }
    }
}
