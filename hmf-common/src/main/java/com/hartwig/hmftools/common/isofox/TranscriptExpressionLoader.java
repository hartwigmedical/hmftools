package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.RnaCommon;
import com.hartwig.hmftools.common.rna.TranscriptExpressionFile;

public final class TranscriptExpressionLoader
{
    public static Map<String,Double> loadTranscriptExpression(final String sampleDir, final String sampleId) throws IOException
    {
        Map<String,Double> transcriptTpms = Maps.newHashMap();

        String filename = TranscriptExpressionFile.generateFilename(sampleDir, sampleId);

        List<String> lines = Files.readAllLines(Paths.get(filename));

        String fileDelim = inferFileDelimiter(filename);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);

        int transNameIndex = fieldsIndexMap.get(FLD_TRANS_NAME);
        int tpmIndex = fieldsIndexMap.get(FLD_ADJ_TPM);

        for(String line : lines.subList(1, lines.size()))
        {
            final String[] items = line.split(fileDelim, -1);

            final String transName = items[transNameIndex];
            double tpm = Double.parseDouble(items[tpmIndex]);

            transcriptTpms.put(transName, tpm);
        }

        return transcriptTpms;
    }
}
