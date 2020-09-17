package com.hartwig.hmftools.common.fusion;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableDisruption
{
    public abstract int svId();

    @NotNull
    public abstract String chromosome();

    public abstract int orientation();

    public abstract int strand();

    @NotNull
    public abstract String chrBand();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double junctionCopyNumber();

    public abstract int exonUp();

    public abstract int exonDown();

    public abstract double undisruptedCopyNumber();

    public static List<ReportableDisruption> from(final List<LinxBreakend> linxBreakends)
    {
        return linxBreakends.stream()
                .filter(x -> x.reportedDisruption())
                .map(x -> ImmutableReportableDisruption.builder()
                        .svId(x.svId())
                        .chromosome(x.chromosome())
                        .orientation(x.orientation())
                        .strand(x.strand())
                        .chrBand(x.chrBand())
                        .gene(x.gene())
                        .type(x.type())
                        .junctionCopyNumber(x.junctionCopyNumber())
                        .exonUp(x.exonUp())
                        .exonDown(x.exonDown())
                        .undisruptedCopyNumber(x.undisruptedCopyNumber())
                        .build())
                .collect(Collectors.toList());
    }

}
