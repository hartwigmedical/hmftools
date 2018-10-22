package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.disruption.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    public static final FieldBuilder<?> LOCATION_FIELD = field("location", String.class);
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
        return new FieldBuilder<?>[] { LOCATION_FIELD, GENE_FIELD, RANGE_FIELD, TYPE_FIELD, COPIES_FIELD, GENE_MIN_COPIES,
                GENE_MAX_COPIES };
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull FittedPurityStatus fitStatus,
            @NotNull List<ReportableGeneDisruption> disruptions) {
        final DRDataSource dataSource = new DRDataSource(LOCATION_FIELD.getName(),
                GENE_FIELD.getName(),
                RANGE_FIELD.getName(),
                TYPE_FIELD.getName(),
                COPIES_FIELD.getName(),
                GENE_MIN_COPIES.getName(),
                GENE_MAX_COPIES.getName());

        for (ReportableGeneDisruption disruption : sort(disruptions)) {
            dataSource.add(disruption.location(),
                    disruption.gene(),
                    disruption.range(),
                    disruption.type().name(),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, ploidyToCopiesString(disruption.ploidy())),
                    Strings.EMPTY,
                    Strings.EMPTY);
            // TODO (KODU): Propagate copy number info here - DEV-548
            //                    PatientReportFormat.correctValueForFitStatus(fitStatus, String.valueOf(disruption.geneMinCopies())),
            //                    PatientReportFormat.correctValueForFitStatus(fitStatus, String.valueOf(disruption.geneMaxCopies())));
        }

        return dataSource;
    }

    @NotNull
    private static List<ReportableGeneDisruption> sort(@NotNull List<ReportableGeneDisruption> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String locationAndGene1 = zeroPrefixed(disruption1.location()) + disruption1.gene();
            String locationAndGene2 = zeroPrefixed(disruption2.location()) + disruption2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return disruption1.firstAffectedExon() - disruption2.firstAffectedExon();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String zeroPrefixed(@NotNull String location) {
        // KODU: First remove q or p arm if present.
        int armStart = location.indexOf("q");
        if (armStart < 0) {
            armStart = location.indexOf("p");
        }

        String chromosome = armStart > 0 ? location.substring(0, armStart) : location;

        try {
            int chromosomeIndex = Integer.valueOf(chromosome);
            if (chromosomeIndex < 10) {
                return "0" + location;
            } else {
                return location;
            }
        } catch (NumberFormatException exception) {
            return location;
        }
    }
}
