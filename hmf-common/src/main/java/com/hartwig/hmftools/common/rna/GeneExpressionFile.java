package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public final class GeneExpressionFile
{
    public static final String GENE_EXPRESSION_FILE_ID = "gene_data.tsv";

    public static final String FLD_SPLICED_FRAGS = "SplicedFragments";
    public static final String FLD_UNSPLICED_FRAGS = "UnsplicedFragments";
    public static final String FLD_ADJ_TPM = "AdjTPM";
    public static final String FLD_MEDIAN_TPM_CANCER = "MedianTpmCancer";
    public static final String FLD_PERC_TPM_CANCER = "PercentileCancer";
    public static final String FLD_MEDIAN_TPM_COHORT = "MedianTpmCohort";
    public static final String FLD_PERC_TPM_COHORT = "PercentileCohort";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + GENE_EXPRESSION_FILE_ID;
    }

    public static List<GeneExpression> read(final String filename)
    {
        try
        {
            List<GeneExpression> geneExpressions = Lists.newArrayList();

            List<String> lines = Files.readAllLines(Paths.get(filename));

            String fileDelim = inferFileDelimiter(filename);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);

            int geneIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int tpmIndex = fieldsIndexMap.get(FLD_ADJ_TPM);
            int splicedFragsIndex = fieldsIndexMap.get(FLD_SPLICED_FRAGS);
            int unsplicedFragsIndex = fieldsIndexMap.get(FLD_UNSPLICED_FRAGS);
            Integer medCancerIndex = fieldsIndexMap.get(FLD_MEDIAN_TPM_CANCER);
            Integer percCancerIndex = fieldsIndexMap.get(FLD_PERC_TPM_CANCER);
            Integer medCohortIndex = fieldsIndexMap.get(FLD_MEDIAN_TPM_COHORT);
            Integer percCohortIndex = fieldsIndexMap.get(FLD_PERC_TPM_COHORT);

            for(String line : lines.subList(1, lines.size()))
            {
                String[] values = line.split(fileDelim, -1);

                geneExpressions.add(ImmutableGeneExpression.builder()
                        .geneName(values[geneIndex])
                        .tpm(Double.parseDouble(values[tpmIndex]))
                        .splicedFragments(Integer.parseInt(values[splicedFragsIndex]))
                        .unsplicedFragments(Integer.parseInt(values[unsplicedFragsIndex]))
                        .medianTpmCancer(getCohortValue(values, medCancerIndex))
                        .percentileCancer(getCohortValue(values, percCancerIndex))
                        .medianTpmCohort(getCohortValue(values, medCohortIndex))
                        .percentileCohort(getCohortValue(values, percCohortIndex))
                        .build());
            }

            return geneExpressions;
        }
        catch(IOException e)
        {
            RNA_LOGGER.error("failed to load Isofox gene data file({}): {}", filename, e.toString());
            return null;
        }
    }

    private static double getCohortValue(final String[] values, Integer fieldIndex)
    {
        if(fieldIndex == null)
            return 0;

        return Double.parseDouble(values[fieldIndex]);
    }
}
