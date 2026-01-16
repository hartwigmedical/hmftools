package com.hartwig.hmftools.datamodel.finding;

import java.io.IOException;
import java.nio.file.Path;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConversionUtil
{
    public static void orangeJsonToFindingsJson(@NotNull Path findingsJson, @NotNull Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        FindingRecord findingRecord = FindingRecordFactory.fromOrangeJsonWithTranscriptFile(orangeJson, clinicalTranscriptsTsv, driverGeneTsv);
        FindingsJson.getInstance().write(findingRecord, findingsJson);
    }
}
