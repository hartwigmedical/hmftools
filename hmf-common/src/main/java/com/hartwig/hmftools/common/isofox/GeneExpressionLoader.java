package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

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
    private static final String FLD_TPM = "AdjTPM";

    private static final String FLD_SPLICED_FRAGS = "SplicedFragments";
    private static final String FLD_UNSPLICED_FRAGS = "UnsplicedFragments";

    @NotNull
    public static List<GeneExpression> loadGeneExpression(
            final String isofoxGeneDataCsv,
            final GeneExpressionDistributionData cohortData,
            final  String cancerType) throws IOException
    {
        List<GeneExpression> geneExpressions = Lists.newArrayList();

        List<String> lines = Files.readAllLines(Paths.get(isofoxGeneDataCsv));

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), RnaCommon.DELIMITER);

        for(String line : lines.subList(1, lines.size()))
        {
            final String[] items = line.split(RnaCommon.DELIMITER, -1);

            final String geneId = items[fieldsIndexMap.get(FLD_GENE_ID)];

            double tpm = Double.parseDouble(items[fieldsIndexMap.get(FLD_TPM)]);

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
