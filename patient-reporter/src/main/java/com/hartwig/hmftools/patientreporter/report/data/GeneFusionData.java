package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.ploidyToCopiesString;

import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneFusionData {

    public abstract String geneStart();

    public abstract String geneContextStart();

    public abstract String geneEnd();

    public abstract String geneContextEnd();

    public abstract String copies();

    @NotNull
    public static GeneFusionData from(@NotNull final GeneFusion fusion) {
        final Transcript upstream = fusion.upstreamLinkedAnnotation();
        final Transcript downstream = fusion.downstreamLinkedAnnotation();

        return ImmutableGeneFusionData.builder()
                .geneStart(upstream.geneName())
                .geneContextStart(exonDescription(upstream, true))
                .geneEnd(downstream.geneName())
                .geneContextEnd(exonDescription(downstream, false))
                .copies(ploidyToCopiesString(fusionPloidy(fusion)))
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