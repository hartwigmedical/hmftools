package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;

import org.jetbrains.annotations.NotNull;

public final class CheckPurpleQuality {

    private CheckPurpleQuality() {
    }

public static boolean checkHasReliablePurity(@NotNull PurityContext purityContext) {
        return purityContext.method() != FittedPurityMethod.NO_TUMOR;
    }

    public static boolean checkHasReliableQuality(@NotNull PurpleQC purpleQC) {
        return purpleQC.status().contains(PurpleQCStatus.PASS);
    }
}
