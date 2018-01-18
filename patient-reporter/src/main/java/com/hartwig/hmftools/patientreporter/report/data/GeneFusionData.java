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
                .type(upstream.variant().type().toString())
                .start(positionString(upstream.geneAnnotation()))
                .geneStart(upstream.geneName())
                .geneContextStart(exonDescription(upstream, true))
                .transcriptStart(upstream.transcriptId())
                .end(positionString(downstream.geneAnnotation()))
                .geneEnd(downstream.geneName())
                .geneContextEnd(exonDescription(downstream, false))
                .transcriptEnd(downstream.transcriptId())
                .vaf(formatNullablePercent(alleleFrequency(upstream.geneAnnotation())) + " " + formatNullablePercent(alleleFrequency(
                        downstream.geneAnnotation())))
                .cosmicURL(fusion.cosmicURL())
                .build();

    }
}