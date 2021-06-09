package com.hartwig.hmftools.protect.purple;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.purple.cnchromosome.ChromosomeArmKey;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleData {

    double purity();

    boolean hasReliablePurity();

    boolean hasReliableQuality();

    double ploidy();

    double microsatelliteIndelsPerMb();

    @NotNull
    MicrosatelliteStatus microsatelliteStatus();

    double tumorMutationalBurdenPerMb();

    int tumorMutationalLoad();

    @NotNull
    TumorMutationalStatus tumorMutationalLoadStatus();

    @NotNull
    List<ReportableVariant> reportableSomaticVariants();

    @NotNull
    List<SomaticVariant> unreportedSomaticExonicVariants();

    @NotNull
    List<ReportableVariant> reportableGermlineVariants();

    @NotNull
    List<ReportableGainLoss> reportableGainsLosses();

    @NotNull
    List<ReportableGainLoss> unreportedGainsLosses();

    @NotNull
    Map<ChromosomeArmKey, Double> cnPerChromosome();
}
