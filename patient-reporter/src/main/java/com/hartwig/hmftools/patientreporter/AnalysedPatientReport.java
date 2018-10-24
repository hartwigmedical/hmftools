package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.chordclassifier.ChordAnalysis;
import com.hartwig.hmftools.patientreporter.disruption.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;

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

    @NotNull
    public abstract FittedPurityStatus fitStatus();

    public abstract double impliedPurity();

    @NotNull
    public abstract List<EvidenceItem> evidenceItems();

    @NotNull
    public abstract List<EnrichedSomaticVariant> somaticVariants();

    @NotNull
    public abstract List<DriverCatalog> somaticVariantDriverCatalog();

    public abstract double microsatelliteIndelsPerMb();

    public abstract List<ChordAnalysis> chordValue();

    public abstract int tumorMutationalLoad();

    public abstract double tumorMutationalBurden();

    public abstract boolean hasGermlineAnalysis();

    @NotNull
    public abstract List<GermlineVariant> germlineVariants();

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
    public abstract String logoRVA();
}
