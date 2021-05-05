package com.hartwig.hmftools.common.metrics;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public final class WGSMetricsFile {

    private static final String MEAN_COVERAGE_COLUMN = "MEAN_COVERAGE";
    private static final String SD_COVERAGE_COLUMN = "SD_COVERAGE";
    private static final String MEDIAN_COVERAGE_COLUMN = "MEDIAN_COVERAGE";
    private static final String MAD_COVERAGE_COLUMN = "MAD_COVERAGE";
    private static final String PCT_EXC_ADAPTER_COLUMN = "PCT_EXC_ADAPTER";
    private static final String PCT_EXC_MAPQ_COLUMN = "PCT_EXC_MAPQ";
    private static final String PCT_EXC_DUPE_COLUMN = "PCT_EXC_DUPE";
    private static final String PCT_EXC_UNPAIRED_COLUMN = "PCT_EXC_UNPAIRED";
    private static final String PCT_EXC_BASEQ_COLUMN = "PCT_EXC_BASEQ";
    private static final String PCT_EXC_OVERLAP_COLUMN = "PCT_EXC_OVERLAP";
    private static final String PCT_EXC_CAPPED_COLUMN = "PCT_EXC_CAPPED";
    private static final String PCT_EXC_TOTAL_COLUMN = "PCT_EXC_TOTAL";

    private static final String COVERAGE_1X_COLUMN = "PCT_1X";
    private static final String COVERAGE_10X_COLUMN = "PCT_10X";
    private static final String COVERAGE_20X_COLUMN = "PCT_20X";
    private static final String COVERAGE_30X_COLUMN = "PCT_30X";
    private static final String COVERAGE_60X_COLUMN = "PCT_60X";

    private WGSMetricsFile() {
    }

    @NotNull
    public static WGSMetrics read(@NotNull String metricsPath) throws IOException {
        WGSMetricsLines lines = WGSMetricsLines.fromFile(metricsPath);

        // These 2 columns do not exist in older versions
        String pctExcAdapter = lines.findValueByHeader(PCT_EXC_ADAPTER_COLUMN);
        String coverage1x = lines.findValueByHeader(COVERAGE_1X_COLUMN);

        return ImmutableWGSMetrics.builder()
                .meanCoverage(Double.parseDouble(lines.findValueByHeader(MEAN_COVERAGE_COLUMN)))
                .sdCoverage(Double.parseDouble(lines.findValueByHeader(SD_COVERAGE_COLUMN)))
                .medianCoverage(Integer.parseInt(lines.findValueByHeader(MEDIAN_COVERAGE_COLUMN)))
                .madCoverage(Integer.parseInt(lines.findValueByHeader(MAD_COVERAGE_COLUMN)))
                .pctExcAdapter(pctExcAdapter != null ? Double.parseDouble(pctExcAdapter) : null)
                .pctExcMapQ(Double.parseDouble(lines.findValueByHeader(PCT_EXC_MAPQ_COLUMN)))
                .pctExcDupe(Double.parseDouble(lines.findValueByHeader(PCT_EXC_DUPE_COLUMN)))
                .pctExcUnpaired(Double.parseDouble(lines.findValueByHeader(PCT_EXC_UNPAIRED_COLUMN)))
                .pctExcBaseQ(Double.parseDouble(lines.findValueByHeader(PCT_EXC_BASEQ_COLUMN)))
                .pctExcOverlap(Double.parseDouble(lines.findValueByHeader(PCT_EXC_OVERLAP_COLUMN)))
                .pctExcCapped(Double.parseDouble(lines.findValueByHeader(PCT_EXC_CAPPED_COLUMN)))
                .pctExcTotal(Double.parseDouble(lines.findValueByHeader(PCT_EXC_TOTAL_COLUMN)))
                .coverage1xPercentage(coverage1x != null ? Double.parseDouble(coverage1x) : null)
                .coverage10xPercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_10X_COLUMN)))
                .coverage20xPercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_20X_COLUMN)))
                .coverage30xPercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_30X_COLUMN)))
                .coverage60xPercentage(Double.parseDouble(lines.findValueByHeader(COVERAGE_60X_COLUMN)))
                .build();
    }
}
