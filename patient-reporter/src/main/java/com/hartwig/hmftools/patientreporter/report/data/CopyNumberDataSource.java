package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class CopyNumberDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    public static final FieldBuilder<?> BAND_FIELD = field("band", String.class);
    public static final FieldBuilder<?> COPY_NUMBER_TYPE_FIELD = field("copynumber_type", String.class);
    public static final FieldBuilder<?> COPY_NUMBER_FIELD = field("copynumber", String.class);

    private CopyNumberDataSource() {
    }

    @NotNull
    public static JRDataSource fromCopyNumbers(@NotNull final List<CopyNumberReport> copyNumbers) {
        final DRDataSource copyNumberDatasource =
                new DRDataSource(CHROMOSOME_FIELD.getName(), BAND_FIELD.getName(), GENE_FIELD.getName(), COPY_NUMBER_TYPE_FIELD.getName(),
                        COPY_NUMBER_FIELD.getName());

        for (final CopyNumberReport copyNumber : copyNumbers) {
            copyNumberDatasource.add(copyNumber.chromosome(), copyNumber.chromosomeBand(), copyNumber.gene(), copyNumber.description(),
                    Integer.toString(copyNumber.copyNumber()));
        }
        return copyNumberDatasource;
    }

    @NotNull
    public static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, BAND_FIELD, GENE_FIELD, COPY_NUMBER_TYPE_FIELD, COPY_NUMBER_FIELD };
    }
}
