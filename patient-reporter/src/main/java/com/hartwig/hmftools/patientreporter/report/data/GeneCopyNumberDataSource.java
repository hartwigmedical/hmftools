package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneCopyNumberDataSource {

    public static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> DESCRIPTION_FIELD = field("description", String.class);
    public static final FieldBuilder<?> COPY_NUMBER_FIELD = field("copynumber", String.class);

    private GeneCopyNumberDataSource() {
    }

    @NotNull
    public static JRDataSource fromCopyNumbers(@NotNull final List<GeneCopyNumber> copyNumbers) {
        final DRDataSource copyNumberDatasource =
                new DRDataSource(POSITION_FIELD.getName(), GENE_FIELD.getName(), DESCRIPTION_FIELD.getName(),
                        COPY_NUMBER_FIELD.getName());

        for (final GeneCopyNumber copyNumber : copyNumbers) {
            copyNumberDatasource.add(copyNumber.chromosomalPosition(), copyNumber.gene(), copyNumber.alteration().name(),
                    Integer.toString((int) Math.round(copyNumber.minCopyNumber())));
        }
        return copyNumberDatasource;
    }

    @NotNull
    public static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { POSITION_FIELD, GENE_FIELD, DESCRIPTION_FIELD, COPY_NUMBER_FIELD };
    }
}
