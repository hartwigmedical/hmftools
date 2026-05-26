package com.hartwig.hmftools.finding;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.util.CandidateToReportableTransformer;
import com.hartwig.hmftools.finding.util.CopyNumberRoundingTransformer;
import com.hartwig.hmftools.finding.util.ErrorTransformer;
import com.hartwig.hmftools.finding.util.LowPurityTransformer;
import com.hartwig.hmftools.finding.util.NoGermlineTransformer;
import com.hartwig.hmftools.finding.util.PTOTransformer;
import com.hartwig.hmftools.finding.util.ReportedOnlyTransformer;
import com.hartwig.hmftools.finding.util.TransformUtil;

import org.jetbrains.annotations.Nullable;

@SuppressWarnings("unused")
public class ConversionUtil
{
    public static FindingRecord orangeRecordToFindingRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, @Nullable Gender gender) throws IOException
    {
        return convert(FindingRecordFactory.fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv, gender));
    }

    private static FindingRecord convert(FindingRecord findingRecord)
    {
        return TransformUtil.listTransformer(List.of(ErrorTransformer::transform,
                LowPurityTransformer::transform,
                PTOTransformer::transform,
                NoGermlineTransformer::transform,
                CandidateToReportableTransformer::transform,
                CopyNumberRoundingTransformer::transform,
                ReportedOnlyTransformer::transform
        )).apply(findingRecord);
    }
}
