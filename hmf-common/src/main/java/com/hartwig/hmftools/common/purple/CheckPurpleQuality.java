package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.purity.PurityContext;

import org.jetbrains.annotations.NotNull;

public final class CheckPurpleQuality {

    private CheckPurpleQuality() {
    }

    public static boolean checkHasReliablePurity(@NotNull PurityContext purityContext) {
        return purityContext.method() != FittedPurityMethod.NO_TUMOR;
    }
}
