package com.hartwig.hmftools.patientreporter.panel;

import com.hartwig.hmftools.patientreporter.PatientReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PanelFailReport implements PatientReport {
}
