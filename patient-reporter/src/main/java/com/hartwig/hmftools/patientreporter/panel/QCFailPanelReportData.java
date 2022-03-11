package com.hartwig.hmftools.patientreporter.panel;

import com.hartwig.hmftools.patientreporter.ReportData;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class QCFailPanelReportData implements ReportData {
}