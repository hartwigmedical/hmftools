package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;

import java.io.File;

public final class GeneExpressionFile
{
    public static final String GENE_EXPRESSION_FILE_ID = "gene_data.csv";

    public static final String FLD_SPLICED_FRAGS = "SplicedFragments";
    public static final String FLD_UNSPLICED_FRAGS = "UnsplicedFragments";
    public static final String FLD_TPM = "AdjTPM";

    public static final String TRANSCRIPT_EXPRESSION_FILE_ID = "transcript_data.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + ISF_FILE_ID + GENE_EXPRESSION_FILE_ID;
    }
}
