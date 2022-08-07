package com.hartwig.hmftools.common.linx;

import java.util.List;

import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxData
{

    @NotNull
    List<LinxSvAnnotation> allStructuralVariants();

    @NotNull
    List<LinxFusion> allFusions();

    @NotNull
    List<LinxFusion> reportableFusions();

    @NotNull
    List<LinxBreakend> allBreakends();

    @NotNull
    List<GeneDisruption> reportableGeneDisruptions();

    @NotNull
    List<HomozygousDisruption> homozygousDisruptions();

    @NotNull
    List<LinxDriver> drivers();

    @NotNull
    List<LinxGermlineSv> allGermlineDisruptions();

    @NotNull
    List<LinxGermlineSv> reportableGermlineDisruptions();
}
