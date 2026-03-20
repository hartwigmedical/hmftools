package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingsJson;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil
{
    public static void orangeJsonToFindingsJson(Path findingsJson, Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, Gender gender) throws IOException
    {
        FindingRecord
                findingRecord = FindingRecordFactory.fromOrangeJsonWithTranscriptFile(orangeJson, clinicalTranscriptsTsv, driverGeneTsv, gender);
        new FindingsJson().write(findingRecord, findingsJson);
    }
}
