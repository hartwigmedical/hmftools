package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SequencedPatientReport implements PatientReport {
    @NotNull
    @Override
    public abstract SampleReport sampleReport();

    @NotNull
    public abstract List<VariantReport> variants();

    @NotNull
    public abstract List<GeneFusionData> geneFusions();

    @NotNull
    public abstract List<GeneDisruptionData> geneDisruptions();

    @NotNull
    public abstract List<CopyNumberReport> copyNumbers();

    public abstract int mutationalLoad();

    @NotNull
    abstract String purplePurity();

    @NotNull
    public String impliedPurityString() {
        return purplePurity();
    }

    @NotNull
    public abstract List<Alteration> civicAlterations();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @NotNull
    @Override
    public abstract String signaturePath();
}
