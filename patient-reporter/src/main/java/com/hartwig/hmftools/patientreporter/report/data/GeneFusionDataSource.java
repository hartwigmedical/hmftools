package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneFusionDataSource {

    public static final FieldBuilder<?> GENES_FIELD = field("genes", String.class);
    public static final FieldBuilder<?> GENE_CONTEXT = field("gene context", String.class);
    public static final FieldBuilder<?> GENE_TRANSCRIPT_FIELD = field("gene_transcript", String.class);
    public static final FieldBuilder<?> PARTNER_CONTEXT_FIELD = field("partner_context", String.class);
    public static final FieldBuilder<?> PARTNER_TRANSCRIPT_FIELD = field("partner_transcript", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    private static final FieldBuilder<?> COSMIC_URL = field("cosmic url", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull List<GeneFusionData> fusions) {
        final DRDataSource dataSource = new DRDataSource(GENES_FIELD.getName(),
                GENE_CONTEXT.getName(),
                GENE_TRANSCRIPT_FIELD.getName(),
                PARTNER_CONTEXT_FIELD.getName(),
                PARTNER_TRANSCRIPT_FIELD.getName(),
                COPIES_FIELD.getName(),
                COSMIC_URL.getName());

        fusions.forEach(fusion -> dataSource.add(fusion.geneStart() + " - " + fusion.geneEnd(),
                fusion.geneContextStart(),
                fusion.geneStartTranscript(),
                fusion.geneContextEnd(),
                fusion.geneEndTranscript(),
                fusion.copies(),
                fusion.cosmicURL()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { GENES_FIELD, GENE_TRANSCRIPT_FIELD, PARTNER_TRANSCRIPT_FIELD, GENE_CONTEXT,/* PARTNER_GENE_FIELD,*/
                PARTNER_CONTEXT_FIELD, COPIES_FIELD, /*URL_TEXT,*/ COSMIC_URL };
    }

    @NotNull
    public static AbstractSimpleExpression<String> cosmicHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return data.getValue(COSMIC_URL.getName());
            }
        };
    }

    @NotNull
    public static AbstractSimpleExpression<String> transcriptUrl(@NotNull final FieldBuilder<?> transcriptField) {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(transcriptField.getName());
            }
        };
    }
}
