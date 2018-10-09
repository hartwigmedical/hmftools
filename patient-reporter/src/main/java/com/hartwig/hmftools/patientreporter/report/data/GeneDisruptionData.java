package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.apache.logging.log4j.util.Strings;
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

    public abstract String affectedRange();

    public abstract String type();

    public abstract String disruptionCopies();

    public abstract String geneMinCopyNumber();

    public abstract String geneMaxCopyNumber();

    @NotNull
    public static GeneDisruptionData from(@NotNull final GeneDisruption disruption) {
        Transcript transcript = disruption.linkedAnnotation();
        GeneAnnotation gene = transcript.parent();
        boolean upstream = gene.variant().orientation(gene.isStart()) > 0;
        String affectedRange = exonDescription(transcript) + (upstream ? " Upstream" : " Downstream");

        return ImmutableGeneDisruptionData.builder()
                .chromosome(gene.variant().chromosome(gene.isStart()))
                .chromosomeBand(gene.karyotypeBand())
                .gene(gene.geneName())
                .affectedRange(affectedRange)
                .type(gene.variant().type().name())
                .disruptionCopies(ploidyToCopiesString(gene.variant().ploidy()))
                .geneMinCopyNumber(Strings.EMPTY)
                .geneMaxCopyNumber(Strings.EMPTY)
                .build();
    }
}
