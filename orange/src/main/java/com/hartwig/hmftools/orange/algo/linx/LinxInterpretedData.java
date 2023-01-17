package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LinxInterpretedData {

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

    @Nullable
    public abstract List<LinxSvAnnotation> allGermlineStructuralVariants();

    @Nullable
    public abstract List<LinxBreakend> allGermlineBreakends();

    @Nullable
    public abstract List<LinxGermlineSv> allGermlineDisruptions();

    @Nullable
    public abstract List<LinxGermlineSv> reportableGermlineDisruptions();
}
