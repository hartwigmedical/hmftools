package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ReportableGeneFusionFactory {

    private ReportableGeneFusionFactory() {
    }

    @NotNull
    public static List<ReportableGeneFusion> fusionConvertToReportable(@NotNull List<Fusion> fusions) {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList();
        for (Fusion fusion : fusions) {
            reportableFusions.add(ImmutableReportableGeneFusion.builder()
                    .geneStart(fusion.geneUp())
                    .geneContextStart(context(fusion.regionTypeUp(), fusion.exonUp(), false))
                    .geneStartTranscript(fusion.transcriptUp())
                    .geneEnd(fusion.geneDown())
                    .geneContextEnd(context(fusion.regionTypeDown(), fusion.exonDown(), true))
                    .geneEndTranscript(fusion.transcriptDown())
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

    private static double fusionPloidy(double downstreamPloidy, double upstreamPloidy) {
        assert Math.abs(upstreamPloidy - downstreamPloidy) < 1E-10;

        return upstreamPloidy;
    }
}
