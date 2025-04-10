package com.hartwig.hmftools.common.sage;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

public final class SageCommon
{
    public static final String SAGE_FILE_ID = ".sage";

    public static String generateBqrFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + SAGE_FILE_ID + ".bqr.tsv";
    }

    public static String generateBqrPlotFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + SAGE_FILE_ID + ".bqr.png";
    }
}
