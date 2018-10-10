package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SomaticVariantAnalysis {

    @NotNull
    public abstract List<EnrichedSomaticVariant> variantsToReport();

    @NotNull
    public abstract List<DriverCatalog> driverCatalog();

    public abstract double microsatelliteIndelsPerMb();

    public abstract int tumorMutationalLoad();

    public abstract double tumorMutationalBurden();

    @NotNull
    public abstract List<EvidenceItem> actionableVariantsReport();

    @NotNull
    public abstract Map<EnrichedSomaticVariant, VariantEvidenceItems> evidence();
}
