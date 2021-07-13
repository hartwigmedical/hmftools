package com.hartwig.hmftools.common.purple;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleData {

    @NotNull
    PurpleQC qc();

    boolean hasReliableQuality();

    @NotNull
    FittedPurityMethod fittedPurityMethod();

    boolean wholeGenomeDuplication();

    boolean hasReliablePurity();

    double purity();

    double minPurity();

    double maxPurity();

    double ploidy();

    double minPloidy();

    double maxPloidy();

    double microsatelliteIndelsPerMb();

    @NotNull
    MicrosatelliteStatus microsatelliteStatus();

    double tumorMutationalBurdenPerMb();

    int tumorMutationalLoad();

    int svTumorMutationalBurden();

    @NotNull
    TumorMutationalStatus tumorMutationalLoadStatus();

    @NotNull
    List<ReportableVariant> reportableSomaticVariants();

    @NotNull
    List<SomaticVariant> unreportedSomaticVariants();

    @NotNull
    List<ReportableVariant> reportableGermlineVariants();

    @NotNull
    List<ReportableGainLoss> reportableGainsLosses();

    @NotNull
    List<ReportableGainLoss> unreportedGainsLosses();

    @NotNull
    List<CnPerChromosomeArmData> cnPerChromosome();
}
