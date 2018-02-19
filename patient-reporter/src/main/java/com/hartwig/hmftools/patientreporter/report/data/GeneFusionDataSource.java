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

    public static final FieldBuilder<?> FUSION_FIELD = field("field", String.class);
    public static final FieldBuilder<?> FIVE_GENE_CONTEXT_FIELD = field("five_gene context", String.class);
    public static final FieldBuilder<?> FIVE_TRANSCRIPT_FIELD = field("five_transcript", String.class);
    public static final FieldBuilder<?> THREE_GENE_CONTEXT_FIELD = field("three_gene_context", String.class);
    public static final FieldBuilder<?> THREE_TRANSCRIPT_FIELD = field("three_transcript", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> COSMIC_URL_TEXT = field("cosmic url text", String.class);
    private static final FieldBuilder<?> COSMIC_URL = field("cosmic url", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull List<GeneFusionData> fusions) {
        final DRDataSource dataSource = new DRDataSource(FUSION_FIELD.getName(),
                FIVE_GENE_CONTEXT_FIELD.getName(),
                FIVE_TRANSCRIPT_FIELD.getName(),
                THREE_GENE_CONTEXT_FIELD.getName(),
                THREE_TRANSCRIPT_FIELD.getName(),
                COPIES_FIELD.getName(),
                COSMIC_URL_TEXT.getName(),
                COSMIC_URL.getName());

        fusions.forEach(fusion -> dataSource.add(fusion.geneStart() + " - " + fusion.geneEnd(),
                fusion.geneContextStart(),
                fusion.geneStartTranscript(),
                fusion.geneContextEnd(),
                fusion.geneEndTranscript(),
                fusion.copies(),
                fusion.cosmicURL().isEmpty() ? "" : "Link",
                fusion.cosmicURL()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { FUSION_FIELD, FIVE_TRANSCRIPT_FIELD, THREE_TRANSCRIPT_FIELD, FIVE_GENE_CONTEXT_FIELD,
                THREE_GENE_CONTEXT_FIELD, COPIES_FIELD, COSMIC_URL_TEXT, COSMIC_URL };
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
