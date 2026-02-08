package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;

public final class TranscriptExpressionFile
{
    public static final String TRANSCRIPT_EXPRESSION_FILE_ID = "transcript_data.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + TRANSCRIPT_EXPRESSION_FILE_ID;
    }
}
