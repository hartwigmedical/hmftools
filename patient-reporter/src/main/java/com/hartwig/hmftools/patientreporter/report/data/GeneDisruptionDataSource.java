package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    public static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    public static final FieldBuilder<?> CHROMOSOME_BAND_FIELD = field("chromosome_band", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> AFFECTED_RANGE_FIELD = field("affected_range", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> GENE_MIN_COPIES = field("gene_min_copies", String.class);
    public static final FieldBuilder<?> GENE_MAX_COPIES = field("gene_max_copies", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull FittedPurityStatus fitStatus, @NotNull List<GeneDisruption> disruptions) {
        final DRDataSource dataSource = new DRDataSource(CHROMOSOME_FIELD.getName(),
                CHROMOSOME_BAND_FIELD.getName(),
                GENE_FIELD.getName(),
                AFFECTED_RANGE_FIELD.getName(),
                TYPE_FIELD.getName(),
                COPIES_FIELD.getName(),
                GENE_MIN_COPIES.getName(),
                GENE_MAX_COPIES.getName());

        final List<GeneDisruptionData> disruptionData =
                disruptions.stream().sorted(disruptionComparator()).map(GeneDisruptionData::from).collect(Collectors.toList());

        disruptionData.forEach(disruption -> dataSource.add(disruption.chromosome(),
                disruption.chromosomeBand(),
                disruption.gene(),
                disruption.affectedRange(),
                disruption.type(),
                PatientReportFormat.correctValueForFitStatus(fitStatus, disruption.disruptionCopies()),
                disruption.geneMinCopyNumber(),
                disruption.geneMaxCopyNumber()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, CHROMOSOME_BAND_FIELD, GENE_FIELD, AFFECTED_RANGE_FIELD, TYPE_FIELD, COPIES_FIELD,
                GENE_MIN_COPIES, GENE_MAX_COPIES};
    }

    @NotNull
    private static Comparator<GeneDisruption> disruptionComparator() {
        return (disruption1, disruption2) -> {
            String location1 = chromosomalLocation(disruption1);
            String location2 = chromosomalLocation(disruption2);

            if (location1.equals(location2)) {
                return disruption1.linkedAnnotation().exonUpstream() - disruption2.linkedAnnotation().exonUpstream();
            } else {
                return location1.compareTo(location2);
            }
        };
    }

    @NotNull
    private static String chromosomalLocation(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = disruption.linkedAnnotation().parent();
        String chromosome = zeroPrefixed(gene.variant().chromosome(gene.isStart()));
        return chromosome + ":" + disruption.linkedAnnotation().parent().karyotypeBand() + ":" + gene.geneName();
    }

    @NotNull
    private static String zeroPrefixed(@NotNull String chromosome) {
        try {
            int chromosomeIndex = Integer.valueOf(chromosome);
            if (chromosomeIndex < 10) {
                return "0" + chromosome;
            } else {
                return chromosome;
            }
        } catch (NumberFormatException exception) {
            return chromosome;
        }
    }
}
