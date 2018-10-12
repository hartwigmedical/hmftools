package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    public static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> RANGE_FIELD = field("range", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> GENE_MIN_COPIES = field("gene_min_copies", String.class);
    public static final FieldBuilder<?> GENE_MAX_COPIES = field("gene_max_copies", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, GENE_FIELD, RANGE_FIELD, TYPE_FIELD, COPIES_FIELD, GENE_MIN_COPIES,
                GENE_MAX_COPIES };
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull FittedPurityStatus fitStatus, @NotNull List<GeneDisruption> disruptions) {
        final DRDataSource dataSource = new DRDataSource(CHROMOSOME_FIELD.getName(),
                GENE_FIELD.getName(),
                RANGE_FIELD.getName(),
                TYPE_FIELD.getName(),
                COPIES_FIELD.getName(),
                GENE_MIN_COPIES.getName(),
                GENE_MAX_COPIES.getName());

        for (GeneDisruption disruption : sort(disruptions)) {
            dataSource.add(chromosomeField(disruption),
                    gene(disruption).geneName(),
                    rangeField(disruption),
                    typeField(disruption),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, copiesField(disruption)),
                    Strings.EMPTY,
                    Strings.EMPTY);
        }

        return dataSource;
    }

    @NotNull
    private static List<GeneDisruption> sort(@NotNull List<GeneDisruption> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String location1 = chromosomalLocationKey(disruption1);
            String location2 = chromosomalLocationKey(disruption2);

            if (location1.equals(location2)) {
                return disruption1.linkedAnnotation().exonUpstream() - disruption2.linkedAnnotation().exonUpstream();
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String chromosomalLocationKey(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        return zeroPrefixed(gene.variant().chromosome(gene.isStart())) + gene.karyotypeBand() + gene.geneName();
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

    @NotNull
    private static String chromosomeField(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        return gene.variant().chromosome(gene.isStart()) + gene.karyotypeBand();
    }

    @NotNull
    private static String rangeField(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        boolean upstream = gene.variant().orientation(gene.isStart()) > 0;

        return exonDescription(disruption.linkedAnnotation()) + (upstream ? " Upstream" : " Downstream");
    }

    @NotNull
    private static String typeField(@NotNull GeneDisruption disruption) {
        return gene(disruption).variant().type().name();
    }

    @NotNull
    private static String copiesField(@NotNull GeneDisruption disruption) {
        return ploidyToCopiesString(gene(disruption).variant().ploidy());
    }

    @NotNull
    private static GeneAnnotation gene(@NotNull GeneDisruption disruption) {
        return disruption.linkedAnnotation().parent();
    }
}
