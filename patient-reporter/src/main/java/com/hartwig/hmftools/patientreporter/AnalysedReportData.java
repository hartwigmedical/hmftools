package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnalysedReportData implements ReportData {

    @NotNull
    public abstract ActionabilityAnalyzer actionabilityAnalyzer();

    @NotNull
    public abstract GermlineReportingModel germlineReportingModel();

    @NotNull
    public abstract SummaryModel summaryModel();
}
