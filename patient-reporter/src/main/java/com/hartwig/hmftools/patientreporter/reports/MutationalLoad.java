package com.hartwig.hmftools.patientreporter.reports;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class MutationalLoad {

    private MutationalLoad() {
    }

    public static int calculate(@NotNull List<SomaticVariant> variants) {
        return 0;
    }
}
