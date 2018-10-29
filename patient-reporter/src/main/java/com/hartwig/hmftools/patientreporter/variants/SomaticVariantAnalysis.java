package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SomaticVariantAnalysis {

    @NotNull
    public abstract List<ReportableSomaticVariant> reportableSomaticVariants();

    @NotNull
    public abstract Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant();

    public abstract double microsatelliteIndelsPerMb();

    public abstract int tumorMutationalLoad();

    public abstract double tumorMutationalBurden();
}
