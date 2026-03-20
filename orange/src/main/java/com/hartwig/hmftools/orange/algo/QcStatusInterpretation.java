package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.datamodel.isofox.RnaQCStatus.FAIL_LOW_COVERAGE;
import static com.hartwig.hmftools.datamodel.purple.PurpleQCStatus.FAIL_CONTAMINATION;
import static com.hartwig.hmftools.datamodel.purple.PurpleQCStatus.FAIL_TINC;

import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;

public final class QcStatusInterpretation
{
    private static final List<PurpleQCStatus> FAIL_STATUS_CODES = List.of(
            PurpleQCStatus.FAIL_NO_TUMOR, FAIL_CONTAMINATION, PurpleQCStatus.FAIL_TINC);

    public static boolean hasPurpleFail(final PurpleQC purpleQC)
    {
        return FAIL_STATUS_CODES.stream().anyMatch(x -> purpleQC.status().contains(x));
    }

    public static boolean hasRnaFail(final IsofoxRecord isofoxRecord)
    {
        return isofoxRecord.summary().qcStatus().stream().anyMatch(x -> x == FAIL_LOW_COVERAGE);
    }

    public static boolean hasTumorContaminated(final PurpleQC purpleQC)
    {
        return purpleQC.status().contains(FAIL_CONTAMINATION) || purpleQC.status().contains(FAIL_TINC);
    }

    public static boolean isFailNoTumor(final PurpleQC purpleQC)
    {
        return purpleQC.status().contains(PurpleQCStatus.FAIL_NO_TUMOR);
    }
}
