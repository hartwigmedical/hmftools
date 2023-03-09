package com.hartwig.hmftools.datamodel.linx;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LinxRecord {

    @NotNull
    public abstract List<LinxSvAnnotation> allSomaticStructuralVariants();

    @NotNull
    public abstract List<LinxFusion> allSomaticFusions();

    @NotNull
    public abstract List<LinxFusion> reportableSomaticFusions();

    @NotNull
    public abstract List<LinxFusion> additionalSuspectSomaticFusions();

    @NotNull
    public abstract List<LinxBreakend> allSomaticBreakends();

    @NotNull
    public abstract List<LinxBreakend> reportableSomaticBreakends();

    @NotNull
    public abstract List<LinxBreakend> additionalSuspectSomaticBreakends();

    @NotNull
    public abstract List<HomozygousDisruption> somaticHomozygousDisruptions();

}
