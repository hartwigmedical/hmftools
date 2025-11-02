package com.hartwig.hmftools.datamodel.linx;

import java.util.List;
import java.util.Optional;

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
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxSvAnnotation> allSomaticStructuralVariants();

    @Nullable
    List<LinxSvAnnotation> allGermlineStructuralVariants();

    @NotNull
    List<LinxFusion> allSomaticFusions();

    @Gson.Ignore
    @NotNull
    default List<LinxFusion> reportableSomaticFusions()
    {
        return allSomaticFusions().stream().filter(LinxFusion::reported).toList();
    }

    @NotNull
    List<LinxBreakend> allSomaticBreakends();

    @Nullable
    List<LinxBreakend> allGermlineBreakends();

    @Gson.Ignore
    @NotNull
    default List<LinxBreakend> reportableSomaticBreakends()
    {
        return allSomaticBreakends().stream().filter(LinxBreakend::reported).toList();
    }

    @Gson.Ignore
    @Nullable
    default List<LinxBreakend> reportableGermlineBreakends()
    {
        return Optional.ofNullable(allGermlineBreakends())
                .map(o -> o.stream().filter(LinxBreakend::reported).toList())
                .orElse(null);
    }

    @NotNull
    List<LinxHomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxHomozygousDisruption> germlineHomozygousDisruptions();
}
