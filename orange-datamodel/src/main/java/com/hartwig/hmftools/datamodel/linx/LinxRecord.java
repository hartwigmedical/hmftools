package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LinxRecord {

    @NotNull
    public abstract List<LinxSvAnnotation> allSomaticStructuralVariants();

    @Nullable
    public abstract List<LinxSvAnnotation> allGermlineStructuralVariants();

    @NotNull
    public abstract List<LinxFusion> allSomaticFusions();

    @NotNull
    public abstract List<LinxFusion> reportableSomaticFusions();

    @NotNull
    public abstract List<LinxFusion> additionalSuspectSomaticFusions();

    @NotNull
    public abstract List<LinxBreakend> allSomaticBreakends();

    @Nullable
    public abstract List<LinxBreakend> allGermlineBreakends();

    @NotNull
    public abstract List<LinxBreakend> reportableSomaticBreakends();

    @Nullable
    public abstract List<LinxBreakend> reportableGermlineBreakends();

    @NotNull
    public abstract List<LinxBreakend> additionalSuspectSomaticBreakends();

    @NotNull
    public abstract List<HomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    public abstract List<HomozygousDisruption> germlineHomozygousDisruptions();
}
