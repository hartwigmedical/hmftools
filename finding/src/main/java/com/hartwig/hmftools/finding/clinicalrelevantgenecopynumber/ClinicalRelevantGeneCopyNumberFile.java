package com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.jetbrains.annotations.NotNull;

public class ClinicalRelevantGeneCopyNumberFile
{

    private ClinicalRelevantGeneCopyNumberFile()
    {
    }

    private enum Column
    {
        Gene
    }

    @NotNull
    public static ClinicalRelevantGeneCopyNumberModel buildFromTsv(@NotNull Path path) throws IOException
    {

        try(DelimFileReader reader = new DelimFileReader(Files.newBufferedReader(path)))
        {
            List<String> clinicalRelevantGeneCopyNumber = reader.stream()
                    .map(row -> row.get(ClinicalRelevantGeneCopyNumberFile.Column.Gene))
                    .toList();
            return new ClinicalRelevantGeneCopyNumberModel(clinicalRelevantGeneCopyNumber);
        }
    }
}