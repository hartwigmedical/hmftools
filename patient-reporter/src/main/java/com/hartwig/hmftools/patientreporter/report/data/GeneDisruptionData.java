package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.util.PatientReportFormat.positionString;

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

    public abstract String geneContext();

    public abstract String copies();

    @NotNull
    public static GeneDisruptionData from(@NotNull final GeneDisruption disruption) {
        final Transcript transcript = disruption.linkedAnnotation();
        final GeneAnnotation gene = transcript.parent();
        final int variantOrientation = gene.variant().orientation(gene.isStart());

        return ImmutableGeneDisruptionData.builder()
                .position(positionString(gene))
                .gene(disruption.linkedAnnotation().geneName())
                .geneContext(exonDescription(transcript, true) + (variantOrientation > 0 ? " Upstream" : " Downstream"))
                .copies("1")
                .build();
    }

}
