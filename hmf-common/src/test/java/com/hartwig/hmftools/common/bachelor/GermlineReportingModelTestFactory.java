package com.hartwig.hmftools.common.bachelor;

import java.util.Map;

import com.hartwig.hmftools.common.bachelor.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

public final class GermlineReportingModelTestFactory {

    private GermlineReportingModelTestFactory() {
    }

    @NotNull
    public static GermlineReportingModel buildFromMap(@NotNull Map<String, Boolean> germlineGenesAndNotifyMap) {
        return new GermlineReportingModel(germlineGenesAndNotifyMap);
    }
}
