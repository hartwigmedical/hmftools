package com.hartwig.hmftools.common.purple;

import java.util.List;

import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
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

    boolean hasReliablePurity();

    double purity();

    double minPurity();

    double maxPurity();

    double ploidy();

    double minPloidy();

    double maxPloidy();

    boolean wholeGenomeDuplication();

    double microsatelliteIndelsPerMb();

    @NotNull
    MicrosatelliteStatus microsatelliteStatus();

    double tumorMutationalBurdenPerMb();

    int tumorMutationalLoad();

    @NotNull
    TumorMutationalStatus tumorMutationalLoadStatus();

    int svTumorMutationalBurden();

    @NotNull
    List<ReportableVariant> reportableSomaticVariants();

    @NotNull
    List<SomaticVariant> unreportedSomaticVariants();

    @NotNull
    List<ReportableVariant> reportableGermlineVariants();

    @NotNull
    List<SomaticVariant> unreportedGermlineVariants();

    @NotNull
    List<ReportableGainLoss> reportableSomaticGainsLosses();

    @NotNull
    List<ReportableGainLoss> unreportedSomaticGainsLosses();

    @NotNull
    List<GermlineDeletion> reportableGermlineDeletions();

    @NotNull
    List<GermlineDeletion> unreportedGermlineDeletions();

    @NotNull
    List<CnPerChromosomeArmData> cnPerChromosome();
}
