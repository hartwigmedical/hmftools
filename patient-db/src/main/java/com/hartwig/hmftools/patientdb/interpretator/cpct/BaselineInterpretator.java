package com.hartwig.hmftools.patientdb.interpretator.cpct;

import org.jetbrains.annotations.Nullable;

public class BaselineInterpretator {

    private BaselineInterpretator() {

    }

    @Nullable
    public static String interpretationCarcinomaPrimaryTumorLocation(@Nullable String primaryTumorLocationCarcinoma,
            @Nullable String primaryTumorLocationOther) {
        // We always read additional info if there is any, see also DEV-1713
        if (primaryTumorLocationCarcinoma != null && primaryTumorLocationOther != null && !primaryTumorLocationOther.isEmpty()) {
            primaryTumorLocationCarcinoma = primaryTumorLocationCarcinoma + " + " + primaryTumorLocationOther;
        }
        return primaryTumorLocationCarcinoma;
    }

    @Nullable
    public static String setFinalPrimaryTumorLocation(@Nullable String  primaryTumorLocationCarcinoma, @Nullable String primaryTumorLocationSelcrit, boolean useCarcinomaForm) {
        // We prefer carcinoma form over sel crit form. See also DEV-540
        return useCarcinomaForm ? primaryTumorLocationCarcinoma : primaryTumorLocationSelcrit;
    }
}
