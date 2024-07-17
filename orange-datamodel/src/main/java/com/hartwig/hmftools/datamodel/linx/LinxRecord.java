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
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxSvAnnotation> allSomaticStructuralVariants();

    @Nullable
    List<LinxSvAnnotation> allGermlineStructuralVariants();

    @NotNull
    List<LinxFusion> allSomaticFusions();

    @NotNull
    List<LinxFusion> reportableSomaticFusions();

    @NotNull
    List<LinxFusion> additionalSuspectSomaticFusions();

    @NotNull
    List<LinxFusion> additionalViableSomaticFusions();

    @NotNull
    List<LinxBreakend> allSomaticBreakends();

    @Nullable
    List<LinxBreakend> allGermlineBreakends();

    @NotNull
    List<LinxBreakend> reportableSomaticBreakends();

    @Nullable
    List<LinxBreakend> reportableGermlineBreakends();

    @NotNull
    List<LinxBreakend> additionalSuspectSomaticBreakends();

    @NotNull
    List<LinxHomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxHomozygousDisruption> germlineHomozygousDisruptions();
}
