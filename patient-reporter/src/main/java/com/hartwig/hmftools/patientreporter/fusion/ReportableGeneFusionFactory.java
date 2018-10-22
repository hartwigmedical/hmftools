package com.hartwig.hmftools.patientreporter.fusion;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableGeneFusionFactory {

    private ReportableGeneFusionFactory() {
    }

    @NotNull
    public static List<ReportableGeneFusion> toReportableGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList();
        for (GeneFusion fusion : fusions) {
            final Transcript upstream = fusion.upstreamLinkedAnnotation();
            final Transcript downstream = fusion.downstreamLinkedAnnotation();

            reportableFusions.add(ImmutableReportableGeneFusion.builder()
                    .geneStart(upstream.geneName())
                    .geneContextStart(exonDescription(upstream))
                    .geneStartTranscript(upstream.transcriptId())
                    .geneEnd(downstream.geneName())
                    .geneContextEnd(exonDescription(downstream))
                    .geneEndTranscript(downstream.transcriptId())
                    .copies(ploidyToCopiesString(fusionPloidy(fusion)))
                    .source(fusion.primarySource())
                    .build());
        }

        return reportableFusions;
    }

    @Nullable
    private static Double fusionPloidy(@NotNull GeneFusion fusion) {
        Double upstreamPloidy = fusion.upstreamLinkedAnnotation().parent().variant().ploidy();
        Double downstreamPloidy = fusion.downstreamLinkedAnnotation().parent().variant().ploidy();

        if (upstreamPloidy == null || downstreamPloidy == null) {
            return null;
        }

        assert upstreamPloidy.equals(downstreamPloidy);

        return upstreamPloidy;
    }
}
