package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class ReportableGeneFusionFactory {

    private ReportableGeneFusionFactory() {
    }

    @NotNull
    public static List<ReportableGeneFusion> convert(@NotNull List<Fusion> fusions) {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList();
        for (Fusion fusion : fusions) {
            reportableFusions.add(ImmutableReportableGeneFusion.builder()
                    .geneStart(fusion.geneUp())
                    .geneContextStart(context(fusion.regionTypeUp(), fusion.exonUp(), false))
                    .geneTranscriptStart(fusion.transcriptUp())
                    .geneEnd(fusion.geneDown())
                    .geneContextEnd(context(fusion.regionTypeDown(), fusion.exonDown(), true))
                    .geneTranscriptEnd(fusion.transcriptDown())
                    .ploidy(fusionPloidy(fusion.ploidyDown(), fusion.ploidyUp()))
                    .source(fusion.primarySource())
                    .build());
        }
        return reportableFusions;
    }

    @NotNull
    private static String context(@NotNull String regionType, int exon, boolean isEnd) {
        switch (regionType) {
            case "Upstream":
                return "Promoter Region";
            case "Exonic":
                return String.format("Exon %d", exon);
            case "Intronic":
                return String.format("Intron %d", isEnd ? exon - 1 : exon);
        }

        return String.format("ERROR: %s", regionType);
    }

    @Nullable
    private static Double fusionPloidy(@Nullable Double downstreamPloidy, @Nullable Double upstreamPloidy) {
        if (downstreamPloidy == null && upstreamPloidy == null) {
            return null;
        }

        assert downstreamPloidy != null && upstreamPloidy != null;
        assert Math.abs(upstreamPloidy - downstreamPloidy) < 1E-10;

        return upstreamPloidy;
    }
}
