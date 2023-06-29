package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;

public final class GeneExpressionFile
{
    public static final String GENE_EXPRESSION_FILE_ID = "gene_data.csv";

    public static final String FLD_SPLICED_FRAGS = "SplicedFragments";
    public static final String FLD_UNSPLICED_FRAGS = "UnsplicedFragments";
    public static final String FLD_ADJ_TPM = "AdjTPM";

    public static final String TRANSCRIPT_EXPRESSION_FILE_ID = "transcript_data.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + GENE_EXPRESSION_FILE_ID;
    }
}
