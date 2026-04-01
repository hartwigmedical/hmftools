package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingsJson;
import com.hartwig.hmftools.finding.util.ErrorConverter;
import com.hartwig.hmftools.finding.util.FindingRecordConverterUtil;
import com.hartwig.hmftools.finding.util.LowPurityConverter;
import com.hartwig.hmftools.finding.util.NoGermlineConverter;
import com.hartwig.hmftools.finding.util.PTOConverter;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil
{
    public static void orangeJsonToFindingsJson(Path findingsJson, Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, Gender gender) throws IOException
    {
        FindingRecord
                findingRecord =
                FindingRecordFactory.fromOrangeJsonWithTranscriptFile(orangeJson, clinicalTranscriptsTsv, driverGeneTsv, gender);
        findingRecord =
                FindingRecordConverterUtil.listConverter(List.of(ErrorConverter::convert,
                                LowPurityConverter::convert,
                                PTOConverter::convert,
                                NoGermlineConverter::convert))
                        .apply(findingRecord);
        new FindingsJson().write(findingRecord, findingsJson);
    }
}
