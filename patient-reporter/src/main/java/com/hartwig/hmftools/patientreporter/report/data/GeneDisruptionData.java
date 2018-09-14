package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

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

    public abstract String chromosome();

    public abstract String chromosomeBand();

    public abstract String gene();

    public abstract String geneContext();

    public abstract String type();

    public abstract String copies();

    @NotNull
    public static GeneDisruptionData from(@NotNull final GeneDisruption disruption) {
        final Transcript transcript = disruption.linkedAnnotation();
        final GeneAnnotation gene = transcript.parent();
        // TODO (KODU): Add upstream/downstream annotation
        //        final boolean upstream = gene.variant().orientation(gene.isStart()) > 0;
        final String geneContext = exonDescription(transcript); //+ (upstream ? " Upstream" : " Downstream");

        return ImmutableGeneDisruptionData.builder()
                .chromosome(gene.variant().chromosome(gene.isStart()))
                .gene(gene.geneName())
                .geneContext(geneContext)
                .type(gene.variant().type().name())
                .copies(ploidyToCopiesString(gene.variant().ploidy()))
                .chromosomeBand(gene.karyotypeBand())
                .build();
    }
}
