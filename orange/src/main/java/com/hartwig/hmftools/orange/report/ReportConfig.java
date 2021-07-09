package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ReportConfig {

    boolean reportGermline();

    @Nullable
    EvidenceLevel maxReportingLevel();

}
