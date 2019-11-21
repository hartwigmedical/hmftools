package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;

import org.jetbrains.annotations.NotNull;

public final class CheckPurpleQuality {

    private CheckPurpleQuality() {

    }

    public static boolean checkHasReliablePurity(@NotNull PurityContext purityContext) {
        return purityContext.status() != FittedPurityStatus.NO_TUMOR;
    }

    public static boolean checkHasReliableQuality(@NotNull PurpleQC purpleQC) {
        return purpleQC.status() == PurpleQCStatus.PASS;
    }

}
