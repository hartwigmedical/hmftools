package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.alleleFrequency;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.formatNullablePercent;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.positionString;

import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneFusionData {

    public abstract String type();

    public abstract String start();

    public abstract String geneStart();

    public abstract String geneContextStart();

    public abstract String transcriptStart();

    public abstract String end();

    public abstract String geneEnd();

    public abstract String geneContextEnd();

    public abstract String transcriptEnd();

    public abstract String vaf();

    @Nullable
    public abstract String cosmicURL();

    public static GeneFusionData from(@NotNull final GeneFusion fusion) {
        final Transcript upstream = fusion.upstreamLinkedAnnotation();
        final Transcript downstream = fusion.downstreamLinkedAnnotation();

        return ImmutableGeneFusionData.builder()
                .type(upstream.getVariant().type().toString())
                .start(positionString(upstream.getGeneAnnotation()))
                .geneStart(upstream.getGeneName())
                .geneContextStart(exonDescription(upstream, true))
                .transcriptStart(upstream.getTranscriptId())
                .end(positionString(downstream.getGeneAnnotation()))
                .geneEnd(downstream.getGeneName())
                .geneContextEnd(exonDescription(downstream, false))
                .transcriptEnd(downstream.getTranscriptId())
                .vaf(formatNullablePercent(alleleFrequency(upstream.getGeneAnnotation())) + " " + formatNullablePercent(
                        alleleFrequency(downstream.getGeneAnnotation())))
                .cosmicURL(fusion.cosmicURL())
                .build();

    }
}