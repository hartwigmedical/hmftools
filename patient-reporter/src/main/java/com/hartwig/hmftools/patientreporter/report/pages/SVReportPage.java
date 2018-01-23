package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.baseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.TextFieldBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class SVReportPage {

    abstract SequencedPatientReport report();

    @NotNull
    private TextFieldBuilder<String> title() {
        return cmp.text(Commons.TITLE + " - Structural Variant Report").setStyle(sectionHeaderStyle());
    }

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() throws IOException {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                title(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                supplementDisclaimerSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneFusionTable(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneDisruptionReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                svExplanation(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneFusionExplanation(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneDisruptionExplanation());
    }

    @NotNull
    private static ComponentBuilder<?, ?> supplementDisclaimerSection() {
        final List<String> lines = Lists.newArrayList("This supplement is a prototype for new types of reporting that "
                        + "may or may not eventually end up in the actual sequencing report",
                "Findings should be considered as indicative and need confirmation by other means "
                        + "before being used in clinical decision making.");

        return toList("Disclaimer", lines);
    }

    @NotNull
    private static ComponentBuilder<?, ?> svExplanation() {
        return toList("Details on structural variants",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "Reported variants are only indicative and have NOT been verified via RNA sequencing."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionExplanation() {
        return toList("Details on reported gene fusions",
                Lists.newArrayList("Only intronic in-frame fusions or whole exon deletions are reported.",
                        "The canonical, or otherwise longest transcript validly fused is reported.",
                        "Fusions are restricted to those in the Fusion Gene list curated by COSMIC.",
                        "We additionally select fusions where one partner occurs in the 5' or 3' position in COSMIC >3 times.",
                        "Whole exon deletions are also restricted to this list.",
                        "See http://cancer.sanger.ac.uk/cosmic/fusion for more information."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionExplanation() {
        return toList("Details on reported gene disruptions",
                Lists.newArrayList("Only the canonical transcript of disrupted genes are reported.",
                        "Reported gene disruptions are restricted to those that occur in the HMF Panel."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionReport(@NotNull final SequencedPatientReport report) {
        final int fontSize = 6;
        final ComponentBuilder<?, ?> table;
        if (report.geneDisruptions().size() > 0) {
            table = cmp.subreport(baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                    .fields(GeneDisruptionDataSource.geneDisruptionFields())
                    .columns(col.column("Gene", GeneDisruptionDataSource.GENE_FIELD).setFixedWidth(50),
                            col.column("Transcript", GeneDisruptionDataSource.TRANSCRIPT_FIELD)
                                    .setHyperLink(hyperLink(fieldTranscriptLink(GeneDisruptionDataSource.TRANSCRIPT_FIELD)))
                                    .setStyle(linkStyle().setFontSize(fontSize)),
                            col.column("Position", GeneDisruptionDataSource.POSITION_FIELD),
                            col.column("Gene Context", GeneDisruptionDataSource.SV_GENE_CONTEXT),
                            col.column("Orientation", GeneDisruptionDataSource.SV_ORIENTATION_FIELD),
                            col.column("Partner", GeneDisruptionDataSource.SV_PARTNER_POSITION_FIELD),
                            col.column("Type", GeneDisruptionDataSource.SV_TYPE_FIELD).setFixedWidth(30),
                            col.column("VAF", GeneDisruptionDataSource.SV_VAF).setFixedWidth(30))
                    .setDataSource(GeneDisruptionDataSource.fromGeneDisruptions(report.geneDisruptions())));
        } else {
            table = cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));
        }

        return cmp.verticalList(cmp.text("Gene Disruptions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionTable(@NotNull final SequencedPatientReport report) {
        final int fontSize = 6;
        final ComponentBuilder<?, ?> table;
        if (report.geneFusions().size() > 0) {
            table = cmp.subreport(baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                    .fields(GeneFusionDataSource.geneFusionFields())
                    .columns(col.column("5' Gene", GeneFusionDataSource.GENE_FIELD).setFixedWidth(50),
                            col.column("5' Transcript", GeneFusionDataSource.TRANSCRIPT_FIELD)
                                    .setHyperLink(hyperLink(fieldTranscriptLink(GeneFusionDataSource.TRANSCRIPT_FIELD)))
                                    .setStyle(linkStyle().setFontSize(fontSize)),
                            col.column("5' Position", GeneFusionDataSource.POSITION_FIELD),
                            col.column("5' Gene Context", GeneFusionDataSource.SV_GENE_CONTEXT),
                            col.column("3' Gene", GeneFusionDataSource.SV_PARTNER_GENE_FIELD).setFixedWidth(50),
                            col.column("3' Transcript", GeneFusionDataSource.SV_PARTNER_TRANSCRIPT_FIELD)
                                    .setHyperLink(hyperLink(fieldTranscriptLink(GeneFusionDataSource.SV_PARTNER_TRANSCRIPT_FIELD)))
                                    .setStyle(linkStyle().setFontSize(fontSize)),
                            col.column("3' Position", GeneFusionDataSource.SV_PARTNER_POSITION_FIELD),
                            col.column("3' Gene Context", GeneFusionDataSource.SV_PARTNER_CONTEXT_FIELD),
                            col.column("SV Type", GeneFusionDataSource.SV_TYPE_FIELD).setFixedWidth(30),
                            col.column("VAF", GeneFusionDataSource.SV_VAF).setFixedWidth(30))
                    .setDataSource(GeneFusionDataSource.fromGeneFusions(report.geneFusions())));
        } else {
            table = cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));
        }

        return cmp.verticalList(cmp.text("Gene Fusions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static AbstractSimpleExpression<String> fieldTranscriptLink(@NotNull final FieldBuilder<?> field) {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(field.getName());
            }
        };
    }
}
