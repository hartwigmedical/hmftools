package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingsJson;
import com.hartwig.hmftools.finding.util.FindingRecordConverterUtil;
import com.hartwig.hmftools.finding.util.LowPurityConverter;
import com.hartwig.hmftools.finding.util.PTOConverter;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil
{
    public static void orangeJsonToFindingsJson(Path findingsJson, Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        FindingRecord
                findingRecord =
                FindingRecordFactory.fromOrangeJsonWithTranscriptFile(orangeJson, clinicalTranscriptsTsv, driverGeneTsv);
        findingRecord =
                FindingRecordConverterUtil.listConverter(List.of(LowPurityConverter::convert, PTOConverter::convert)).apply(findingRecord);
        new FindingsJson().write(findingRecord, findingsJson);
    }
}
