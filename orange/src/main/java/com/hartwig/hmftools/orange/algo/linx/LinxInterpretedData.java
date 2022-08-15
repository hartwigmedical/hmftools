package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.linx.GeneDisruption;
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
    public abstract List<LinxSvAnnotation> allStructuralVariants();

    @NotNull
    public abstract List<LinxFusion> allFusions();

    @NotNull
    public abstract List<LinxFusion> reportableFusions();

    @NotNull
    public abstract List<LinxFusion> additionalSuspectFusions();

    @NotNull
    public abstract List<LinxBreakend> allBreakends();

    @NotNull
    public abstract List<GeneDisruption> reportableGeneDisruptions();

    @NotNull
    public abstract List<GeneDisruption> additionalSuspectDisruptions();

    @NotNull
    public abstract List<HomozygousDisruption> homozygousDisruptions();

    @NotNull
    public abstract List<LinxGermlineSv> allGermlineDisruptions();

    @NotNull
    public abstract List<LinxGermlineSv> reportableGermlineDisruptions();
}
