package com.hartwig.hmftools.protect.linx;

import java.util.List;

import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxData {

    @NotNull
    List<LinxFusion> fusions();

    @NotNull
    List<ReportableGeneDisruption> geneDisruptions();

    @NotNull
    List<ReportableHomozygousDisruption> homozygousDisruptions();

    @NotNull
    List<ViralInsertion> viralInsertions();
}
