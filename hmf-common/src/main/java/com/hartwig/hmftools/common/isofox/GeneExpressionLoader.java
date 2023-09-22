package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.common.rna.RnaCommon;

import org.jetbrains.annotations.NotNull;

public final class GeneExpressionLoader
{
    @NotNull
    public static List<GeneExpression> loadGeneExpression(
            final String isofoxGeneDataFile,
            final GeneExpressionDistributionData cohortData,
            final  String cancerType) throws IOException
    {
        List<GeneExpression> geneExpressions = Lists.newArrayList();

        List<String> lines = Files.readAllLines(Paths.get(isofoxGeneDataFile));

        String fileDelim = inferFileDelimiter(isofoxGeneDataFile);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);

        for(String line : lines.subList(1, lines.size()))
        {
            final String[] items = line.split(fileDelim, -1);

            final String geneId = items[fieldsIndexMap.get(FLD_GENE_ID)];

            double tpm = Double.parseDouble(items[fieldsIndexMap.get(FLD_ADJ_TPM)]);

            double medianCancer = cohortData.getTpmMedian(geneId, cancerType);
            double percentileCancer = cohortData.getTpmPercentile(geneId, cancerType, tpm);
            double medianCohort = cohortData.getTpmMedian(geneId, GeneExpressionDistributionData.PAN_CANCER);
            double percentileCohort = cohortData.getTpmPercentile(geneId, GeneExpressionDistributionData.PAN_CANCER, tpm);

            geneExpressions.add(ImmutableGeneExpression.builder()
                    .geneName(items[fieldsIndexMap.get(FLD_GENE_NAME)])
                    .tpm(tpm)
                    .splicedFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_SPLICED_FRAGS)]))
                    .unsplicedFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_UNSPLICED_FRAGS)]))
                    .medianTpmCancer(medianCancer)
                    .percentileCancer(percentileCancer)
                    .medianTpmCohort(medianCohort)
                    .percentileCohort(percentileCohort)
                    .build());
        }

        return geneExpressions;
    }
}
