package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleRecord {

    @NotNull
    public abstract PurpleFit fit();

    @NotNull
    public abstract PurpleCharacteristics characteristics();

    @NotNull
    public abstract List<PurpleDriver> somaticDrivers();

    @Nullable
    public abstract List<PurpleDriver> germlineDrivers();

    @NotNull
    public abstract List<PurpleVariant> allSomaticVariants();

    @NotNull
    public abstract List<PurpleVariant> reportableSomaticVariants();

    @Nullable
    public abstract List<PurpleVariant> allGermlineVariants();

    @Nullable
    public abstract List<PurpleVariant> reportableGermlineVariants();

    @NotNull
    public abstract List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    public abstract List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers();

    @NotNull
    public abstract List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH();

    @NotNull
    public abstract List<PurpleGainLoss> allSomaticGainsLosses();

    @NotNull
    public abstract List<PurpleGainLoss> reportableSomaticGainsLosses();
}
