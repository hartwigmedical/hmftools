package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;


@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityVariantData {

    public abstract String event();

    public abstract String matchingCancerType();

    public abstract String source();

    public abstract String drug();

    public abstract String drugsType();

    public abstract String level();

    public abstract String response();

    @NotNull
    public static ActionabilityVariantData from() {
        return ImmutableActionabilityVariantData.builder()
                .event("a")
                .matchingCancerType("a")
                .source("a")
                .drug("a")
                .drugsType("a")
                .level("a")
                .response("a")
                .build();
    }
}
