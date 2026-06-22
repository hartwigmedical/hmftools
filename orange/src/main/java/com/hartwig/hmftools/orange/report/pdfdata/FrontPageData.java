package com.hartwig.hmftools.orange.report.pdfdata;

import java.util.Map;

import org.jetbrains.annotations.Nullable;

public record FrontPageData(Map<String, String> sampleSummary, @Nullable String qcWarning, Map<String, String> technicalSummary,
                            Map<String, String> driverSummary, Map<String, String> genomeWideFeatures, String circosPlotPath)
{
}
