package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneFusionData {
    @NotNull
    public abstract String geneStart();

    @NotNull
    public abstract String geneContextStart();

    @NotNull
    public abstract String geneStartTranscript();

    @NotNull
    public abstract String geneEnd();

    @NotNull
    public abstract String geneContextEnd();

    @NotNull
    public abstract String geneEndTranscript();

    @NotNull
    public abstract String copies();

    @NotNull
    public abstract String source();

    @NotNull
    public static GeneFusionData from(@NotNull final GeneFusion fusion) {
        final Transcript upstream = fusion.upstreamLinkedAnnotation();
        final Transcript downstream = fusion.downstreamLinkedAnnotation();

        return ImmutableGeneFusionData.builder()
                .geneStart(upstream.geneName())
                .geneContextStart(exonDescription(upstream))
                .geneStartTranscript(upstream.transcriptId())
                .geneEnd(downstream.geneName())
                .geneContextEnd(exonDescription(downstream))
                .geneEndTranscript(downstream.transcriptId())
                .copies(ploidyToCopiesString(fusionPloidy(fusion)))
                .source(fusion.primarySource())
                .build();
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
