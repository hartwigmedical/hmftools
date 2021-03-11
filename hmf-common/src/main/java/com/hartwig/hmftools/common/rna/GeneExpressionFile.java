package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;

import java.io.File;

import org.jetbrains.annotations.NotNull;

public class GeneExpressionFile
{
    public static final String GENE_DATA_FILE_ID = "gene_data.csv";

    public GeneExpressionFile()
    {
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + ISF_FILE_ID + GENE_DATA_FILE_ID;
    }
}
