package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleFit
{
    @NotNull
    PurpleQC qc();

    @NotNull
    PurpleFittedPurityMethod fittedPurityMethod();

    double purity();

    double minPurity();

    double maxPurity();

    double ploidy();

    double minPloidy();

    double maxPloidy();

    default boolean containsTumorCells()
    {
        return !qc().status().contains(PurpleQCStatus.FAIL_NO_TUMOR);
    }
}
