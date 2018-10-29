package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantReport {

    @NotNull
    public abstract List<ReportableGeneFusion> reportableFusions();

    @NotNull
    public abstract List<ReportableGeneDisruption> reportableDisruptions();

    @NotNull
    public abstract Map<GeneFusion, List<EvidenceItem>> evidencePerFusion();
}
