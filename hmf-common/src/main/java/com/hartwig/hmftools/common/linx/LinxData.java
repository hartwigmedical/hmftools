package com.hartwig.hmftools.common.linx;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxData
{
    @NotNull
    List<LinxSvAnnotation> allSomaticStructuralVariants();

    @NotNull
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxFusion> allSomaticFusions();

    @NotNull
    List<LinxFusion> reportableSomaticFusions();

    @NotNull
    List<LinxBreakend> allSomaticBreakends();

    @NotNull
    List<LinxBreakend> reportableSomaticBreakends();

    @NotNull
    List<HomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxSvAnnotation>  allGermlineStructuralVariants();

    @Nullable
    List<LinxBreakend> allGermlineBreakends();

    @Nullable
    List<LinxGermlineSv> allGermlineDisruptions();

    @Nullable
    List<LinxGermlineSv> reportableGermlineDisruptions();
}
