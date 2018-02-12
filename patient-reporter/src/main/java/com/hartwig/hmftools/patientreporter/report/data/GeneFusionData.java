package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.ploidyToCopiesString;

import java.util.List;

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
    public abstract List<Integer> geneStartEntrezIds();

    @NotNull
    public abstract String geneContextStart();

    @NotNull
    public abstract String geneEnd();

    @NotNull
    public abstract List<Integer> geneEndEntrezIds();

    @NotNull
    public abstract String geneContextEnd();

    @NotNull
    public abstract String copies();

    @NotNull
    public abstract String cosmicURL();

    @NotNull
    public static GeneFusionData from(@NotNull final GeneFusion fusion) {
        final Transcript upstream = fusion.upstreamLinkedAnnotation();
        final Transcript downstream = fusion.downstreamLinkedAnnotation();

        return ImmutableGeneFusionData.builder()
                .geneStart(upstream.geneName())
                .geneStartEntrezIds(upstream.parent().entrezIds())
                .geneContextStart(exonDescription(upstream))
                .geneEnd(downstream.geneName())
                .geneEndEntrezIds(downstream.parent().entrezIds())
                .geneContextEnd(exonDescription(downstream))
                .copies(ploidyToCopiesString(fusionPloidy(fusion)))
                .cosmicURL(fusion.cosmicURL())
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