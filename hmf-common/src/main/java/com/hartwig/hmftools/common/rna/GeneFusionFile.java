package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

public final class GeneFusionFile
{
    public static final String UNFILTERED_FUSION_FILE_ID = "fusions.csv";
    public static final String PASS_FUSION_FILE_ID = "pass_fusions.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + PASS_FUSION_FILE_ID;
    }

    public static String generateUnfilteredFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + UNFILTERED_FUSION_FILE_ID;
    }
}
