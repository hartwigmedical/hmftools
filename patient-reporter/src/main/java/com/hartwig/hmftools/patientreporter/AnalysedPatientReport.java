package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnalysedPatientReport implements PatientReport {
    @NotNull
    @Override
    public abstract SampleReport sampleReport();

    public abstract double impliedPurity();

    public abstract boolean hasReliablePurityFit();

    public abstract double averageTumorPloidy();

    @NotNull
    public abstract List<EvidenceItem> tumorSpecificEvidence();

    @NotNull
    public abstract List<ClinicalTrial> clinicalTrials();

    @NotNull
    public abstract List<EvidenceItem> offLabelEvidence();

    @NotNull
    public abstract List<ReportableVariant> somaticVariants();

    public abstract double microsatelliteIndelsPerMb();

    public abstract int tumorMutationalLoad();

    public abstract double tumorMutationalBurden();

    @NotNull
    public abstract ChordAnalysis chordAnalysis();

    @Nullable
    public abstract LimsGermlineReportingChoice germlineOptionPatient();

    public abstract boolean reportableGermlineVariant();

    @NotNull
    public abstract List<GeneCopyNumber> geneCopyNumbers();

    @NotNull
    public abstract List<ReportableGeneFusion> geneFusions();

    @NotNull
    public abstract List<ReportableGeneDisruption> geneDisruptions();

    @NotNull
    public abstract String circosPath();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @NotNull
    @Override
    public abstract String signaturePath();

    @NotNull
    @Override
    public abstract String logoRVAPath();
}
