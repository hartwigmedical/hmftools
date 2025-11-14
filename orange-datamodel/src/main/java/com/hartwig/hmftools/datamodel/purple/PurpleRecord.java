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
    TumorStats tumorStats();

    @NotNull
    PurpleCharacteristics characteristics();

    @NotNull
    List<PurpleDriver> somaticDrivers();

    @Nullable
    List<PurpleDriver> germlineDrivers();

    @NotNull
    List<PurpleVariant> driverSomaticVariants();

    @NotNull
    List<PurpleVariant> otherSomaticVariants();

    @Nullable
    List<PurpleVariant> driverGermlineVariants();

    @Nullable
    List<PurpleVariant> otherGermlineVariants();

    @NotNull
    List<PurpleCopyNumber> somaticCopyNumbers();

    @NotNull
    List<PurpleGeneCopyNumber> somaticGeneCopyNumbers();

    @NotNull
    List<PurpleGainDeletion> driverSomaticGainsDels();

    @Nullable
    List<PurpleGermlineDeletion> otherGermlineDeletions();

    @Nullable
    List<PurpleGainDeletion> driverGermlineDeletions();

    @Nullable
    List<PurpleLossOfHeterozygosity> driverGermlineLossOfHeterozygosities();

    @Nullable
    List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities();
}
