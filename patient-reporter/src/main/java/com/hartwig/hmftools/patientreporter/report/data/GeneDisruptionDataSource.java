package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    public static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> GENE_CONTEXT = field("gene context", String.class);
    public static final FieldBuilder<?> COPIES = field("copies", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull List<GeneDisruptionData> disruptions) {
        final DRDataSource dataSource = new DRDataSource(POSITION_FIELD.getName(),
                GENE_FIELD.getName(),
                TYPE_FIELD.getName(),
                GENE_CONTEXT.getName(),
                COPIES.getName());

        disruptions.forEach(disruption -> dataSource.add(disruption.position(),
                disruption.gene(),
                disruption.type(),
                disruption.geneContext(),
                disruption.copies()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { POSITION_FIELD, GENE_FIELD, TYPE_FIELD, GENE_CONTEXT, COPIES };
    }
}
