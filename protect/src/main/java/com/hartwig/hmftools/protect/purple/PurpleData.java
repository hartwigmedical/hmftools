package com.hartwig.hmftools.protect.purple;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.variants.ReportableVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleData {

    double purity();

    double ploidy();

    double microsatelliteIndelsPerMb();

    @NotNull
    MicrosatelliteStatus microsatelliteStatus();

     double tumorMutationalBurdenPerMb();

    int tumorMutationalLoad();

    @NotNull
    TumorMutationalStatus tumorMutationalLoadStatus();

    @NotNull
    List<GeneCopyNumber> geneCopyNumbers();

    @NotNull
    List<ReportableVariant> somaticVariants();

    @NotNull
    List<ReportableGainLoss> copyNumberAlterations();
}
