package com.hartwig.hmftools.datamodel.purple;

import java.util.List;
import java.util.Optional;

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
    List<PurpleVariant> allSomaticVariants();

    @Gson.Ignore
    @NotNull
    default List<PurpleVariant> reportableSomaticVariants()
    {
        return allSomaticVariants().stream().filter(PurpleVariant::reported).toList();
    }

    @Nullable
    List<PurpleVariant> allGermlineVariants();

    @Gson.Ignore
    @Nullable
    default List<PurpleVariant> reportableGermlineVariants()
    {
        return Optional.ofNullable(allGermlineVariants())
                .map(o -> o.stream().filter(PurpleVariant::reported).toList())
                .orElse(null);
    }

    @NotNull
    List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers();

    @NotNull
    List<PurpleGainDeletion> reportableSomaticGainsDels();

    @Nullable
    List<PurpleGermlineDeletion> allGermlineDeletions();

    @Nullable
    List<PurpleGainDeletion> allGermlineFullDels();

    @Nullable
    List<PurpleGainDeletion> reportableGermlineFullDels();

    @Nullable
    List<PurpleLossOfHeterozygosity> allGermlineLossOfHeterozygosities();

    @Nullable
    List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities();
}
