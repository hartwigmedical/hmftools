package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Fusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableGeneFusionFactory {

    private ReportableGeneFusionFactory() {
    }

    @NotNull
    public static List<ReportableGeneFusion> fusionConvertToReportable(@Nullable List<Fusion> fusions) {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList();
        if (fusions != null) {
            for (Fusion fusion: fusions) {
                reportableFusions.add(ImmutableReportableGeneFusion.builder()
                        .geneStart(fusion.geneUp())
                        .geneContextStart(fusion.biotypeUp())
                        .geneStartTranscript(fusion.transStartUp())
                        .geneEnd(fusion.geneDown())
                        .geneContextEnd(fusion.biotypeUp())
                        .geneStartTranscript(fusion.transEndDown())
                        .ploidy(0)
                        .source(fusion.primarySource())
                        .build());
            }
        }
        return reportableFusions;
    }

    @NotNull
    public static List<ReportableGeneFusion> toReportableGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList();
        for (GeneFusion fusion : fusions) {
            final Transcript upstream = fusion.upstreamTrans();
            final Transcript downstream = fusion.downstreamTrans();

            reportableFusions.add(ImmutableReportableGeneFusion.builder()
                    .geneStart(upstream.geneName())
                    .geneContextStart(exonDescription(upstream))
                    .geneStartTranscript(upstream.StableId)
                    .geneEnd(downstream.geneName())
                    .geneContextEnd(exonDescription(downstream))
                    .geneEndTranscript(downstream.StableId)
                    .ploidy(fusionPloidy(fusion))
                    .source(fusion.primarySource())
                    .build());
        }

        return reportableFusions;
    }

    @NotNull
    private static Double fusionPloidy(@NotNull GeneFusion fusion) {
        Double upstreamPloidy = fusion.upstreamTrans().parent().variant().ploidy();
        Double downstreamPloidy = fusion.downstreamTrans().parent().variant().ploidy();

        if (upstreamPloidy == null || downstreamPloidy == null) {
            // Not sure when ploidy would be null...
            return Double.NaN;
        }

        assert upstreamPloidy.equals(downstreamPloidy);

        return upstreamPloidy;
    }
}
