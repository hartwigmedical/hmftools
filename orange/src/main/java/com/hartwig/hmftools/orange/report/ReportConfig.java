package com.hartwig.hmftools.orange.report;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ReportConfig {

    boolean limitJsonOutput();

    boolean reportGermline();

}
