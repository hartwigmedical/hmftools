package com.hartwig.hmftools.common.linx;

import java.util.List;

import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxData {

    @NotNull
    List<LinxFusion> allFusions();

    @NotNull
    List<LinxFusion> reportableFusions();

    @NotNull
    List<ReportableGeneDisruption> geneDisruptions();

    @NotNull
    List<ReportableHomozygousDisruption> homozygousDisruptions();

    @NotNull
    List<LinxDriver> drivers();

    @NotNull
    List<LinxGermlineSv> allGermlineDisruptions();

    @NotNull
    List<LinxGermlineSv> reportableGermlineDisruptions();
}
