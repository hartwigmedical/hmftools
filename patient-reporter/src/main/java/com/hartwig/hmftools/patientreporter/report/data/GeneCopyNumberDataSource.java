package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneCopyNumberDataSource {

    public static final FieldBuilder<?> CHROMOSOME = field("chromosome", String.class);
    public static final FieldBuilder<?> CHROMOSOME_BAND = field("chromosome band", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> GAIN_OR_LOSS_FIELD = field("gain_or_loss", String.class);
    public static final FieldBuilder<?> COPY_NUMBER_FIELD = field("copynumber", String.class);

    private GeneCopyNumberDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { CHROMOSOME, CHROMOSOME_BAND, GENE_FIELD, GAIN_OR_LOSS_FIELD, COPY_NUMBER_FIELD };
    }

    @NotNull
    public static JRDataSource fromCopyNumbers(@NotNull final List<GeneCopyNumber> copyNumbers, boolean hasReliablePurityFit) {
        final DRDataSource copyNumberDatasource = new DRDataSource(CHROMOSOME.getName(),
                CHROMOSOME_BAND.getName(),
                GENE_FIELD.getName(),
                GAIN_OR_LOSS_FIELD.getName(),
                COPY_NUMBER_FIELD.getName());

        for (GeneCopyNumber copyNumber : sort(copyNumbers)) {
            copyNumberDatasource.add(copyNumber.chromosome(),
                    copyNumber.chromosomeBand(),
                    copyNumber.gene(),
                    type(copyNumber),
                    PatientReportFormat.correctValueForFitReliability(Integer.toString(copyNumber.value()), hasReliablePurityFit));
        }
        return copyNumberDatasource;
    }

    @NotNull
    private static List<GeneCopyNumber> sort(@NotNull List<GeneCopyNumber> geneCopyNumbers) {
        return geneCopyNumbers.stream().sorted((copyNumber1, copyNumber2) -> {
            String location1 = zeroPrefixed(copyNumber1.chromosome() + copyNumber1.chromosomeBand());
            String location2 = zeroPrefixed(copyNumber2.chromosome() + copyNumber2.chromosomeBand());

            if (location1.equals(location2)) {
                return copyNumber1.gene().compareTo(copyNumber2.gene());
            } else {
                return location1.compareTo(location2);
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

    @NotNull
    public static String type(@NotNull GeneCopyNumber geneCopyNumber) {
        if (geneCopyNumber.alteration() == CopyNumberAlteration.GAIN) {
            return "gain";
        } else {
            // KODU: At this point we only have losses and gains.
            assert geneCopyNumber.alteration() == CopyNumberAlteration.LOSS;
            if (geneCopyNumber.maxCopyNumber() < 0.5) {
                return "full loss";
            } else {
                return "partial loss";
            }
        }
    }
}
