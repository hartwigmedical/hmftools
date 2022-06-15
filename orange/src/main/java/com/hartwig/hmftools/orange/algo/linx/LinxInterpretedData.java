package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LinxInterpretedData {

    @NotNull
    public abstract List<LinxFusion> allFusions();

    @NotNull
    public abstract List<LinxFusion> reportableFusions();

    @NotNull
    public abstract List<LinxFusion> additionalSuspectFusions();

    @NotNull
    public abstract List<LinxBreakend> allBreakends();

    @NotNull
    public abstract List<ReportableGeneDisruption> reportableGeneDisruptions();

    @NotNull
    public abstract List<LinxBreakend> additionalSuspectBreakends();

    @NotNull
    public abstract List<ReportableHomozygousDisruption> homozygousDisruptions();

    @NotNull
    public abstract List<LinxGermlineSv> allGermlineDisruptions();

    @NotNull
    public abstract List<LinxGermlineSv> reportableGermlineDisruptions();
}
