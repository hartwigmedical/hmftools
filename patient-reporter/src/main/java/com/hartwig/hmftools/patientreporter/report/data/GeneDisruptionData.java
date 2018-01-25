package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.exonDescription;

import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneDisruptionData {

    public abstract String position();

    public abstract String gene();

    public abstract String type();

    public abstract String geneContext();

    public abstract String copies();

    @NotNull
    public static GeneDisruptionData from(@NotNull final GeneDisruption disruption) {
        final Transcript transcript = disruption.linkedAnnotation();
        final GeneAnnotation gene = transcript.parent();
        final boolean upstream = gene.variant().orientation(gene.isStart()) > 0;

        return ImmutableGeneDisruptionData.builder()
                .position("todo")
                .gene(gene.geneName())
                .type(gene.variant().type().name())
                .geneContext(exonDescription(transcript, upstream) + (upstream ? " Upstream" : " Downstream"))
                .copies("1")
                .build();
    }
}
