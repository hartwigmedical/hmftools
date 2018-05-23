package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;

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
    public static JRDataSource fromCopyNumbers(@NotNull final List<GeneCopyNumber> copyNumbers) {
        final DRDataSource copyNumberDatasource = new DRDataSource(CHROMOSOME.getName(),
                CHROMOSOME_BAND.getName(),
                GENE_FIELD.getName(),
                GAIN_OR_LOSS_FIELD.getName(),
                COPY_NUMBER_FIELD.getName());

        for (final GeneCopyNumber copyNumber : copyNumbers) {
            copyNumberDatasource.add(copyNumber.chromosome(),
                    copyNumber.chromosomeBand(),
                    copyNumber.gene(), type(copyNumber),
                    Integer.toString(copyNumber.value()));
        }
        return copyNumberDatasource;
    }

    @NotNull
    public static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { CHROMOSOME, CHROMOSOME_BAND, GENE_FIELD, GAIN_OR_LOSS_FIELD, COPY_NUMBER_FIELD };
    }

    @NotNull
    private static String type(@NotNull GeneCopyNumber geneCopyNumber) {
        String valueLoss = "";
        if (geneCopyNumber.alteration() == CopyNumberAlteration.GAIN) {
            return "gain";
        } else {
            assert geneCopyNumber.alteration() == CopyNumberAlteration.LOSS;
            if (geneCopyNumber.maxCopyNumber() <= 0.5) {
                valueLoss = "full loss";
            } else if (geneCopyNumber.maxCopyNumber() >= 2.0) {
                valueLoss = "partial loss";
            }
            return valueLoss;
        }
    }
}
