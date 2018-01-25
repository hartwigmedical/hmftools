package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneFusionDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> GENE_CONTEXT = field("gene context", String.class);
    public static final FieldBuilder<?> PARTNER_GENE_FIELD = field("partner_gene", String.class);
    public static final FieldBuilder<?> PARTNER_CONTEXT_FIELD = field("partner_context", String.class);
    public static final FieldBuilder<?> COPIES = field("copies", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull List<GeneFusionData> fusions) {
        final DRDataSource dataSource = new DRDataSource(GENE_FIELD.getName(),
                GENE_CONTEXT.getName(),
                PARTNER_GENE_FIELD.getName(),
                PARTNER_CONTEXT_FIELD.getName(),
                COPIES.getName());

        fusions.forEach(fusion -> dataSource.add(fusion.geneStart(),
                fusion.geneContextStart(),
                fusion.geneEnd(),
                fusion.geneContextEnd(),
                fusion.copies()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, GENE_CONTEXT, PARTNER_GENE_FIELD, PARTNER_CONTEXT_FIELD, COPIES };
    }
}
