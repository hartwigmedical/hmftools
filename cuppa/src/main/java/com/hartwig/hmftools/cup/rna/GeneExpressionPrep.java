package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.RNA;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class GeneExpressionPrep implements CategoryPrep
{
    private final PrepConfig mConfig;
    private static final String FLOAT_FORMAT_LOG_TPM = "%6.3e";

    public GeneExpressionPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return CategoryType.GENE_EXP; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();
        final String filename = mConfig.geneExpressionFile(sampleId);

        try
        {
            String fileDelim = inferFileDelimiter(filename);

            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);
            fileData.remove(0);

            int geneNameCol = fieldsIndexMap.get(FLD_GENE_NAME);
            int adjTPM = fieldsIndexMap.get(FLD_ADJ_TPM);

            for(String line : fileData)
            {
                final String[] items = line.split(fileDelim, -1);

                String geneName = items[geneNameCol];
                double adjTpm = Double.parseDouble(items[adjTPM]);

                double logTpm = log(adjTpm + 1);

                dataItems.add(new DataItem(RNA, ItemType.EXPRESSION, geneName, logTpm, FLOAT_FORMAT_LOG_TPM));
            }

            return dataItems;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to extract category({})", sampleId, categoryType());
            e.printStackTrace();
            return null; // No System.exit(1) allows RNA data to be missing for a sample in multi-sample mode
        }
    }
}
