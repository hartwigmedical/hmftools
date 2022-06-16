package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleInterpretedData {

    @NotNull
    public abstract PurityPloidyFit fit();

    @NotNull
    public abstract PurpleCharacteristics characteristics();

    @NotNull
    public abstract List<SomaticVariant> allSomaticVariants();

    @NotNull
    public abstract List<ReportableVariant> reportableSomaticVariants();

    @NotNull
    public abstract List<ReportableVariant> additionalSuspectSomaticVariants();

    @NotNull
    public abstract List<SomaticVariant> allGermlineVariants();

    @NotNull
    public abstract List<ReportableVariant> reportableGermlineVariants();

    @NotNull
    public abstract List<ReportableVariant> additionalSuspectGermlineVariants();

    @NotNull
    public abstract List<GeneCopyNumber> allSomaticGeneCopyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH();

    @NotNull
    public abstract List<GainLoss> allSomaticGainsLosses();

    @NotNull
    public abstract List<GainLoss> reportableSomaticGainsLosses();

    @NotNull
    public abstract List<GainLoss> nearReportableSomaticGains();

    @NotNull
    public abstract List<GainLoss> additionalSuspectSomaticGainsLosses();

    @NotNull
    public abstract List<GermlineDeletion> allGermlineDeletions();

    @NotNull
    public abstract List<GermlineDeletion> reportableGermlineDeletions();

    @NotNull
    public abstract List<CnPerChromosomeArmData> copyNumberPerChromosome();
}
