package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.StructuralVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientReport {
    @NotNull
    public abstract String sample();

    @NotNull
    public abstract List<VariantReport> variants();

    @NotNull
    public abstract List<StructuralVariantAnalysis.GeneFusion> geneFusions();

    @NotNull
    public abstract List<StructuralVariantAnalysis.GeneDisruption> geneDisruptions();

    @NotNull
    public abstract List<CopyNumberReport> copyNumbers();

    public abstract int mutationalLoad();

    @NotNull
    public abstract String tumorType();

    @Nullable
    abstract Double tumorPercentage();

    @NotNull
    abstract String purplePurity();

    public String tumorPercentageString() {
        return PatientReportFormat.formatNullablePercent(tumorPercentage());
    }

    @NotNull
    public String impliedPurityString() {
        return purplePurity();
    }

    @Nullable
    public abstract String tumorBarcode();

    @Nullable
    public abstract String bloodBarcode();

    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @Nullable
    public abstract LocalDate bloodArrivalDate();
}
