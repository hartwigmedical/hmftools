package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRange;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRangeEvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

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
    public abstract List<EnrichedSomaticVariant> somaticVariants();

    @NotNull
    public abstract List<EvidenceItem> somaticActionabilityVariants();

    @NotNull
    public abstract List<ActionabilityRange> somaticActionabilityVariantsRange();

    @NotNull
    public abstract Map<EnrichedSomaticVariant, VariantEvidenceItems> evidenceItem();

    @NotNull
    public abstract Map<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> evidenceItemRange();

    @NotNull
    public abstract List<DriverCatalog> somaticVariantDriverCatalog();

    public abstract double microsatelliteIndelsPerMb();

    public abstract int tumorMutationalLoad();

    public abstract double tumorMutationalBurden();

    public abstract List<GermlineVariant> germlineVariants();

    @NotNull
    public abstract List<GeneCopyNumber> geneCopyNumbers();

    @NotNull
    public abstract List<GeneFusion> geneFusions();

    @NotNull
    public abstract List<GeneDisruption> geneDisruptions();

    @NotNull
    public abstract String circosPath();

    @Override
    @NotNull
    public abstract Optional<String> comments();

    @NotNull
    @Override
    public abstract String signaturePath();
}
