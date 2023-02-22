package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;

import java.io.File;

public final class TranscriptExpressionFile
{
    public static final String TRANSCRIPT_EXPRESSION_FILE_ID = "transcript_data.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + ISF_FILE_ID + TRANSCRIPT_EXPRESSION_FILE_ID;
    }
}
