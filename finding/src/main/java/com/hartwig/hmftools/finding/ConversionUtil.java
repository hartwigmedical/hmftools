package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.util.CandidateToReportableConverter;
import com.hartwig.hmftools.finding.util.CopyNumberConverter;
import com.hartwig.hmftools.finding.util.CopyNumberRoundingConverter;
import com.hartwig.hmftools.finding.util.ErrorConverter;
import com.hartwig.hmftools.finding.util.FindingRecordConverterUtil;
import com.hartwig.hmftools.finding.util.GainDeletionsFilterConverter;
import com.hartwig.hmftools.finding.util.LowPurityConverter;
import com.hartwig.hmftools.finding.util.NoGermlineConverter;
import com.hartwig.hmftools.finding.util.PTOConverter;
import com.hartwig.hmftools.finding.util.ReportedOnlyConverter;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil {
    public static FindingRecord orangeRecordToFindingRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, @Nullable Gender gender) throws IOException {
        return convert(FindingRecordFactory.fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv, gender));
    }

    private static FindingRecord convert(FindingRecord findingRecord) {
        return FindingRecordConverterUtil.listConverter(List.of(ErrorConverter::convert,
                LowPurityConverter::convert,
                PTOConverter::convert,
                GainDeletionsFilterConverter::convert,
                NoGermlineConverter::convert,
                CandidateToReportableConverter::convert,
                CopyNumberConverter::convert,
                ReportedOnlyConverter::convert
                )).apply(findingRecord);
    }
}
