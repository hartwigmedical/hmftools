package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    public static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    public static final FieldBuilder<?> CHROMOSOME_BAND_FIELD = field("chromosome band", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> GENE_CONTEXT_FIELD = field("gene context", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull FittedPurityStatus fitStatus, @NotNull List<GeneDisruption> disruptions) {
        final DRDataSource dataSource = new DRDataSource(CHROMOSOME_FIELD.getName(),
                CHROMOSOME_BAND_FIELD.getName(),
                GENE_FIELD.getName(),
                GENE_CONTEXT_FIELD.getName(),
                TYPE_FIELD.getName(),
                COPIES_FIELD.getName());

        final List<GeneDisruptionData> disruptionData =
                disruptions.stream().sorted(disruptionComparator()).map(GeneDisruptionData::from).collect(Collectors.toList());

        disruptionData.forEach(disruption -> dataSource.add(disruption.chromosome(),
                disruption.chromosomeBand(),
                disruption.gene(),
                disruption.geneContext(),
                disruption.type(),
                PatientReportFormat.correctValueForFitStatus(fitStatus, disruption.copies())));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, CHROMOSOME_BAND_FIELD, GENE_FIELD, GENE_CONTEXT_FIELD, TYPE_FIELD, COPIES_FIELD };
    }

    @NotNull
    private static Comparator<GeneDisruption> disruptionComparator() {
        return Comparator.comparing(GeneDisruption::linkedAnnotation, Comparator.comparing((Transcript transcript) -> {
            final long startPosition = transcript.parent().variant().start().position();
            // TODO (KODU) : Review if this leads to proper sorting. Idea is that disruptions are sorted by chromosome/band.
            Long codingStart = transcript.codingStart();
            Long codingEnd = transcript.codingEnd();
            if (codingStart != null && codingEnd != null && startPosition >= codingStart && startPosition <= codingEnd) {
                return transcript.parent().variant().start();
            } else {
                return transcript.parent().variant().end();
            }
        }));
    }
}
