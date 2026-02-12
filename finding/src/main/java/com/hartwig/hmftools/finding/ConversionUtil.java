package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.datamodel.finding.FindingRecord;
import com.hartwig.hmftools.datamodel.finding.FindingsJson;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil
{
    public static void orangeJsonToFindingsJson(Path findingsJson, Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        FindingRecord
                findingRecord = FindingRecordFactory.fromOrangeJsonWithTranscriptFile(orangeJson, clinicalTranscriptsTsv, driverGeneTsv);
        new FindingsJson().write(findingRecord, findingsJson);
    }
}
