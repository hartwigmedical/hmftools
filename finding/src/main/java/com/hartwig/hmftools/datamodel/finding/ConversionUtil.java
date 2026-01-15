package com.hartwig.hmftools.datamodel.finding;

import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConversionUtil
{

    public static void orangeJsonToFindingsJson(@NotNull Path findingsJson, @NotNull Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        OrangeRecord orangeRecord = OrangeJson.getInstance().read(orangeJson.toString());
        FindingRecord findingRecord = FindingRecordFactory.fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv);
        FindingsJson.getInstance().write(findingRecord, findingsJson);
    }
}
