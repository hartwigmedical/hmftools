package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.NotNull;

public final class ReportWriterFactory {

    private ReportWriterFactory() {
    }

    @NotNull
    public static ReportWriter createToDiskWriter(@NotNull OrangeConfig config) {
        return new ReportWriter(true, config.outputDir(), config.reportConfig());
    }

    @NotNull
    public static ReportWriter createInMemoryWriter(@NotNull OrangeConfig config) {
        return new ReportWriter(false, null, config.reportConfig());
    }
}
