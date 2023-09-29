package com.hartwig.hmftools.orange.report.interpretation;

import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCInterpretation
{
    public static boolean isFail(@NotNull PurpleQC purpleQC)
    {
        return isFailNoTumor(purpleQC) || isContaminated(purpleQC);
    }

    public static boolean isContaminated(@NotNull PurpleQC purpleQC)
    {
        return purpleQC.status().contains(PurpleQCStatus.FAIL_CONTAMINATION);
    }

    public static boolean isFailNoTumor(@NotNull PurpleQC purpleQC)
    {
        return purpleQC.status().contains(PurpleQCStatus.FAIL_NO_TUMOR);
    }
}
