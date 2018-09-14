package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.CGI;
import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.CIVIC;
import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.COSMIC;
import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.ONCOKB;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneFusionDataSource {

    public static final FieldBuilder<?> FUSION_FIELD = field("fusion", String.class);
    public static final FieldBuilder<?> START_TRANSCRIPT_FIELD = field("five_transcript", String.class);
    public static final FieldBuilder<?> END_TRANSCRIPT_FIELD = field("three_transcript", String.class);
    public static final FieldBuilder<?> START_CONTEXT_FIELD = field("five_gene_context", String.class);
    public static final FieldBuilder<?> END_CONTEXT_FIELD = field("three_gene_context", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> SOURCE_FIELD = field("source", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull FittedPurityStatus fitStatus, @NotNull List<GeneFusion> fusions) {
        final DRDataSource dataSource = new DRDataSource(FUSION_FIELD.getName(),
                START_TRANSCRIPT_FIELD.getName(),
                END_TRANSCRIPT_FIELD.getName(),
                START_CONTEXT_FIELD.getName(),
                END_CONTEXT_FIELD.getName(),
                COPIES_FIELD.getName(),
                SOURCE_FIELD.getName());

        final List<GeneFusionData> fusionsData =
                fusions.stream().sorted(fusionComparator()).map(GeneFusionData::from).collect(Collectors.toList());

        fusionsData.forEach(fusion -> dataSource.add(name(fusion),
                fusion.geneStartTranscript(),
                fusion.geneEndTranscript(),
                fusion.geneContextStart(),
                fusion.geneContextEnd(),
                PatientReportFormat.correctValueForFitStatus(fitStatus, fusion.copies()),
                fusion.source()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { FUSION_FIELD, START_TRANSCRIPT_FIELD, END_TRANSCRIPT_FIELD, START_CONTEXT_FIELD, END_CONTEXT_FIELD,
                COPIES_FIELD, SOURCE_FIELD };
    }

    @NotNull
    private static String name(@NotNull GeneFusionData fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
    }

    @NotNull
    public static AbstractSimpleExpression<String> transcriptUrl(@NotNull final FieldBuilder<?> transcriptField) {
        return new TranscriptExpression(transcriptField);
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                final String source = data.getValue(SOURCE_FIELD.getName());
                switch (source) {
                    case ONCOKB:
                        return "http://oncokb.org/#/";
                    case COSMIC:
                        return "https://cancer.sanger.ac.uk/cosmic";
                    case CGI:
                        return "https://www.cancergenomeinterpreter.org/biomarkers";
                    case CIVIC:
                        return "https://civicdb.org/browse/somaticVariants";
                    default:
                        return "";
                }
            }
        };
    }

    @NotNull
    private static Comparator<GeneFusion> fusionComparator() {
        return Comparator.comparing(fusion -> fusion.upstreamLinkedAnnotation().geneName());
    }
}
