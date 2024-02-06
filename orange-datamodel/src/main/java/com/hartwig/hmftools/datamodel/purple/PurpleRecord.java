package com.hartwig.hmftools.datamodel.purple;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleRecord
{
    @NotNull
    PurpleFit fit();

    @NotNull
    PurpleCharacteristics characteristics();

    @NotNull
    List<PurpleDriver> somaticDrivers();

    @Nullable
    List<PurpleDriver> germlineDrivers();

    @NotNull
    List<PurpleVariant> allSomaticVariants();

    @NotNull
    List<PurpleVariant> reportableSomaticVariants();

    @NotNull
    List<PurpleVariant> additionalSuspectSomaticVariants();

    @Nullable
    List<PurpleVariant> allGermlineVariants();

    @Nullable
    List<PurpleVariant> reportableGermlineVariants();

    @Nullable
    List<PurpleVariant> additionalSuspectGermlineVariants();

    @NotNull
    List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers();

    @NotNull
    List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH();

    @NotNull
    List<PurpleGainLoss> allSomaticGainsLosses();

    @NotNull
    List<PurpleGainLoss> reportableSomaticGainsLosses();

    @NotNull
    List<PurpleGainLoss> nearReportableSomaticGains();

    @NotNull
    List<PurpleGainLoss> additionalSuspectSomaticGainsLosses();

    @Nullable
    List<PurpleGermlineDeletion> allGermlineDeletions();

    @Nullable
    List<PurpleGainLoss> allGermlineFullLosses();

    @Nullable
    List<PurpleGainLoss> reportableGermlineFullLosses();

    @Nullable
    List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities();

    @Nullable
    List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities();
}
