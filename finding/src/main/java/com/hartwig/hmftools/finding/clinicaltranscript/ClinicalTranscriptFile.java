package com.hartwig.hmftools.finding.clinicaltranscript;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTranscriptFile
{
    private enum Column
    {
        Gene,
        Transcript
    }

    private ClinicalTranscriptFile()
    {
    }

    @NotNull
    public static ClinicalTranscriptsModel buildFromTsv(@NotNull OrangeRefGenomeVersion orangeRefGenomeVersion,
            @NotNull Path path) throws IOException
    {
        if(orangeRefGenomeVersion == OrangeRefGenomeVersion.V37)
        {
            try (DelimFileReader reader = new DelimFileReader(Files.newBufferedReader(path)))
            {
                Map<String, String> clinicalTranscriptEntries = reader.stream()
                        .collect(Collectors.toMap(row -> row.get(Column.Gene), row -> row.get(Column.Transcript)));
                return new ClinicalTranscriptsModel(clinicalTranscriptEntries);
            }
        }
        return ClinicalTranscriptsModel.emptyModel();
    }
}