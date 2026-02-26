package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public final class TranscriptExpressionFile
{
    public static final String TRANSCRIPT_EXPRESSION_FILE_ID = "transcript_data.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + TRANSCRIPT_EXPRESSION_FILE_ID;
    }

    public static List<TranscriptExpression> read(String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(Paths.get(filename));

        String fileDelim = inferFileDelimiter(filename);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);

        int transNameIndex = fieldsIndexMap.get(FLD_TRANS_NAME);
        int geneNameIndex = fieldsIndexMap.get(FLD_GENE_NAME);
        int tpmIndex = fieldsIndexMap.get(FLD_ADJ_TPM);

        List<TranscriptExpression> transcriptExpressions = new ArrayList<>();
        for(String line : lines.subList(1, lines.size()))
        {
            final String[] items = line.split(fileDelim, -1);

            final String transName = items[transNameIndex];
            final String geneName = items[geneNameIndex];
            double tpm = Double.parseDouble(items[tpmIndex]);

            transcriptExpressions.add(ImmutableTranscriptExpression.builder()
                    .transcriptName(transName)
                    .geneName(geneName)
                    .tpm(tpm)
                    .build());
        }

        return transcriptExpressions;
    }
}
