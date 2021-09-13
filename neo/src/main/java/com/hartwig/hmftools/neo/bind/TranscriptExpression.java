package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class TranscriptExpression
{
    private final Map<String,Double> mTranscriptExpression;
    private final boolean mValidData;

    public static final String IMMUNE_EXPRESSION_FILE = "trans_exp_file";
    public static final String IMMUNE_EXPRESSION_FILE_CFG = "Immunogenic transcript expression";

    public TranscriptExpression(final String filename)
    {
        mTranscriptExpression = Maps.newHashMap();
        mValidData = loadExpressionData(filename);
    }

    public boolean hasData() { return mValidData; }

    public Double getExpression(final String transName)
    {
        return mTranscriptExpression.get(transName);
    }

    private boolean loadExpressionData(final String filename)
    {
        if(filename == null)
            return false;

        if(!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("transcript expression file({}) not found", filename);
            return false;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            // GeneName,TransName,TPM
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int transcriptIndex = fieldsIndexMap.get("TransName");
            int tpmIndex = fieldsIndexMap.get("TPM");

            for(String line :lines)
            {
                final String[] values = line.split(DELIMITER, -1);

                String transName = values[transcriptIndex];
                double tpm = Double.parseDouble(values[tpmIndex]);

                mTranscriptExpression.put(transName, tpm);
            }

            NE_LOGGER.info("loaded expression for {} transcripts from file({})", mTranscriptExpression.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read transcript expression data file: {}", e.toString());
            return false;
        }

        return true;
    }
}
