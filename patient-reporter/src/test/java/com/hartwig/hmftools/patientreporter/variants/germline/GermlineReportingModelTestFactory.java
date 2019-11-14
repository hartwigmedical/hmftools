package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.Map;


import org.jetbrains.annotations.NotNull;

public final class GermlineReportingModelTestFactory {

    private GermlineReportingModelTestFactory() {
    }

    @NotNull
    public static GermlineReportingModel buildFromMap(@NotNull Map<String, Boolean> germlineGenesAndNotifyMap) {
        return new GermlineReportingModel(germlineGenesAndNotifyMap);
    }
}
