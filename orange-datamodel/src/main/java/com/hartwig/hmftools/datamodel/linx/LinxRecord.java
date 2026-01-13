package com.hartwig.hmftools.datamodel.linx;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxRecord
{
    @NotNull
    List<LinxSvAnnotation> somaticStructuralVariants();

    @Nullable
    List<LinxSvAnnotation> germlineStructuralVariants();

    @NotNull
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxFusion> fusions();

    @NotNull
    List<LinxBreakend> somaticBreakends();

    @Nullable
    List<LinxBreakend> germlineBreakends();

    @NotNull
    List<LinxHomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxHomozygousDisruption> germlineHomozygousDisruptions();
}
