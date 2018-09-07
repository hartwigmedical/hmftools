package com.hartwig.hmftools.patientreporter.variants;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantReport {

    @NotNull
    SomaticVariant variant();

    @NotNull
    String gene();

    int totalReadCount();

    int alleleReadCount();

    @Value.Derived
    default double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }

    @NotNull
    @Value.Derived
    default String readDepthField() {
        return alleleReadCount() + " / " + totalReadCount() + " (" + PatientReportFormat.formatPercent(alleleFrequency()) + ")";
    }

    @NotNull
    String proteinImpact();

    @NotNull
    String proteinImpactType();

    @Nullable
    String knowledgebaseKey();

    @Nullable
    String knowledgebaseUrl();

    @NotNull
    String ploidy();

    double purityAdjustedVAF();

    @NotNull
    @Value.Derived
    default String purityAdjustedVAFField() {
        return PatientReportFormat.formatPercent(purityAdjustedVAF());
    }

    @NotNull
    String clonalityStatus();

    @NotNull
    String wildTypeStatus();

    @NotNull
    String driverStatus();

    @NotNull
    String actionabilityStatus();
}
