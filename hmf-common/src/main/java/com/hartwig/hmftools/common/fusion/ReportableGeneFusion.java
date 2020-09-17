package com.hartwig.hmftools.common.fusion;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGeneFusion {

    @NotNull
    public abstract String geneStart();

    @NotNull
    public abstract String geneContextStart();

    @NotNull
    public abstract String geneTranscriptStart();

    @NotNull
    public abstract String geneEnd();

    @NotNull
    public abstract String geneContextEnd();

    @NotNull
    public abstract String geneTranscriptEnd();

    @Nullable
    public abstract Double junctionCopyNumber();

    public static List<ReportableGeneFusion> from(final List<LinxFusion> linxFusions)
    {
        return linxFusions.stream()
                .filter(x -> x.reported())
                .map(x -> ImmutableReportableGeneFusion.builder()
                        .geneStart(x.geneStart())
                        .geneContextStart(x.geneContextStart())
                        .geneTranscriptStart(x.geneTranscriptStart())
                        .geneEnd(x.geneEnd())
                        .geneContextEnd(x.geneContextEnd())
                        .geneTranscriptEnd(x.geneTranscriptEnd())
                        .junctionCopyNumber(x.junctionCopyNumber())
                        .build())
                .collect(Collectors.toList());
    }
}
