package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneFusionDataSource {

    public static final FieldBuilder<?> FUSION_FIELD = field("field", String.class);
    public static final FieldBuilder<?> START_TRANSCRIPT_FIELD = field("five_transcript", String.class);
    public static final FieldBuilder<?> END_TRANSCRIPT_FIELD = field("three_transcript", String.class);
    public static final FieldBuilder<?> START_CONTEXT_FIELD = field("five_gene context", String.class);
    public static final FieldBuilder<?> END_CONTEXT_FIELD = field("three_gene_context", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> COSMIC_URL_TEXT = field("cosmic_url_text", String.class);
    private static final FieldBuilder<?> COSMIC_URL = field("cosmic_url", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull List<GeneFusionData> fusions) {
        final DRDataSource dataSource = new DRDataSource(FUSION_FIELD.getName(),
                START_TRANSCRIPT_FIELD.getName(),
                END_TRANSCRIPT_FIELD.getName(),
                START_CONTEXT_FIELD.getName(),
                END_CONTEXT_FIELD.getName(),
                COPIES_FIELD.getName(),
                COSMIC_URL_TEXT.getName(),
                COSMIC_URL.getName());

        fusions.forEach(fusion -> dataSource.add(name(fusion),
                fusion.geneStartTranscript(),
                fusion.geneEndTranscript(),
                fusion.geneContextStart(),
                fusion.geneContextEnd(),
                fusion.copies(),
                fusion.cosmicURL().isEmpty() ? Strings.EMPTY : name(fusion),
                fusion.cosmicURL()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { FUSION_FIELD, START_TRANSCRIPT_FIELD, END_TRANSCRIPT_FIELD, START_CONTEXT_FIELD, END_CONTEXT_FIELD,
                COPIES_FIELD, COSMIC_URL_TEXT, COSMIC_URL };
    }

    @NotNull
    private static String name(@NotNull GeneFusionData fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
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
        return new TranscriptExpression(transcriptField);
    }
}
