package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(PDFWriter.class);

    // MIVO: change font to monospace to remove text truncation issue (see gene panel type column for example)
    private static final String FONT = "Times New Roman";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);

    private static final int TEXT_HEADER_INDENT = 30;
    private static final int TEXT_DETAIL_INDENT = 40;
    private static final int LIST_INDENT = 5;
    private static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;
    private static final int DETAIL_TO_DETAIL_VERTICAL_GAP = 4;
    private static final int SECTION_VERTICAL_GAP = 25;
    private static final int PADDING = 1;

    @NotNull
    private final String reportDirectory;
    @NotNull
    @VisibleForTesting
    static final String REPORT_LOGO_PATH = "pdf/hartwig_logo.jpg";

    public PDFWriter(@NotNull final String reportDirectory) {
        this.reportDirectory = reportDirectory;
    }

    @Override
    public void writeSequenceReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException {
        final InputStream logoStream = Resources.asByteSource(Resources.getResource(REPORT_LOGO_PATH)).openStream();
        final JasperReportBuilder reportBuilder = generatePatientReport(report, logoStream, reporterData);
        writeReport(report.sample(), reportBuilder);
        logoStream.close();
    }

    @Override
    public void writeNonSequenceableReport(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage, @NotNull final NotSequenceableReason reason, @NotNull final NotSequenceableStudy study)
            throws IOException, DRException {
        final InputStream logoStream = Resources.asByteSource(Resources.getResource(REPORT_LOGO_PATH)).openStream();
        final JasperReportBuilder reportBuilder =
                generateNotSequenceableReport(sample, tumorType, tumorPercentage, reason, study, logoStream);
        writeReport(sample, reportBuilder);
        logoStream.close();
    }

    private void writeReport(@NotNull final String sample, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        final String fileName = fileName(sample);
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write report as it already exists: " + fileName);
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created patient report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_hmf_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generateNotSequenceableReport(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentageString, @NotNull final NotSequenceableReason reason,
            @NotNull final NotSequenceableStudy study, @NotNull final InputStream logoStream) {
        // @formatter:off
        final ComponentBuilder<?, ?> report =
                cmp.verticalList(
                        mainPageTopSection(sample, tumorType, tumorPercentageString, logoStream),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageNotSequenceableSection(reason, study));
        // @formatter:on

        return report().noData(report);
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report, @NotNull final InputStream logoStream,
            @NotNull final HmfReporterData reporterData) {
        // @formatter:off
        final ComponentBuilder<?, ?> reportMainPage =
                cmp.verticalList(
                        mainPageTopSection(report.sample(), report.tumorType(), report.tumorPercentageString(),
                                logoStream),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageAboutSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantReport(report, reporterData),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberReport(report));

        final ComponentBuilder<?, ?> structuralVariantPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Structural Variant Information")
                                .setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        geneFusionReport(report),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        geneDisruptionReport(report)
                );

        final ComponentBuilder<?, ?> genePanelPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Gene Panel Information")
                                .setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        genePanelSection(reporterData)
                );

        final ComponentBuilder<?, ?> additionalInfoPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Additional Information")
                                .setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantFieldExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        disclaimerSection()
                );

        final ComponentBuilder<?, ?> totalReport =
                cmp.multiPageList()
                        .add(reportMainPage)
                        .newPage()
                        .add(structuralVariantPage)
                        .newPage()
                        .add(genePanelPage)
                        .newPage()
                        .add(additionalInfoPage);
        // @formatter:on

        return report().noData(totalReport);
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageTopSection(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage, @NotNull final InputStream logoStream) {
        // @formatter:off
        final ComponentBuilder<?, ?> mainDiagnosisInfo = cmp.horizontalList(
                cmp.verticalList(
                        cmp.text("Report Date").setStyle(tableHeaderStyle()),
                        cmp.currentDate().setPattern("dd-MMM-yyyy").setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Primary Tumor Location").setStyle(tableHeaderStyle()),
                        cmp.text(tumorType).setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Pathology Tumor Percentage").setStyle(tableHeaderStyle()),
                        cmp.text(tumorPercentage).setStyle(dataTableStyle()))
        );
        return cmp.horizontalList(
                cmp.image(logoStream),
                cmp.verticalList(
                        cmp.text("HMF Sequencing Report - " + sample)
                                .setStyle(fontStyle().bold().setFontSize(14)
                                    .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                                .setHeight(50),
                        mainDiagnosisInfo)
        );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageAboutSection() {
        return toList("About this report",
                Lists.newArrayList("This test is performed for research purpose and is not meant to be used for clinical decision making.",
                        "Additional information on the various fields can be found on the final page of this report.",
                        "For DRUP-specific questions, please contact the DRUP study team at DRUP@nki.nl.",
                        "For other questions, please contact us via info@hartwigmedicalfoundation.nl."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageNotSequenceableSection(@NotNull final NotSequenceableReason reason,
            @NotNull final NotSequenceableStudy study) {
        final String title;
        final String subTitle;
        final String message;

        switch (reason) {
            case LOW_DNA_YIELD: {
                title = "Notification tumor sample on hold for sequencing";
                subTitle = "Insufficient amount of DNA";
                message = "The amount of isolated DNA was <300 ng, which is insufficient for sequencing. "
                        + "This sample is on hold for further processing awaiting optimization of protocols.";
                break;
            }
            case LOW_TUMOR_PERCENTAGE: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Insufficient percentage of tumor cells";
                message = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            default: {
                title = "TITLE";
                subTitle = "SUB_TITLE";
                message = "MESSAGE";
            }
        }

        // @formatter:off
        return cmp.verticalList(
                cmp.text(title).setStyle(tableHeaderStyle().setFontSize(12))
                        .setHeight(20),
                cmp.text(subTitle).setStyle(dataTableStyle().setFontSize(12))
                        .setHeight(20),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(message).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The received biopsies for the tumor sample for this patient were inadequate to obtain a reliable sequencing " +
                                "result. Therefore whole genome sequencing cannot be performed, " +
                                "unless additional fresh tumor material can be provided for a new assessment.")
                        .setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("When possible, please resubmit using the same " + study.studyName() + "-number. In case additional tumor " +
                        "material cannot be provided, please be notified that the patient will not be evaluable " +
                        "for the " + study.studyCode() + " study.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("For questions, please contact us via info@hartwigmedicalfoundation.nl")
                        .setStyle(fontStyle()));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData) {
        final String mutationalLoadAddition =
                "Patients with a mutational load over 140 could be " + "eligible for immunotherapy within DRUP.";

        final String geneMutationAddition = "Marked genes (*) are included in the DRUP study and indicate potential "
                + "eligibility in DRUP. Please note that the marking is NOT based on the specific variant reported for "
                + "this sample, but only on a gene-level.";

        // @formatter:off
        final ComponentBuilder<?, ?> table = report.variants().size() > 0 ?
                cmp.subreport(baseTable().fields(PatientDataSource.variantFields())
                        .columns(
                            col.column("Gene", PatientDataSource.GENE_FIELD).setFixedWidth(50),
                            col.column("Position", PatientDataSource.POSITION_FIELD),
                            col.column("Variant", PatientDataSource.VARIANT_FIELD),
                            col.column("Depth (VAF)", PatientDataSource.DEPTH_VAF_FIELD),
                            col.componentColumn("Predicted Effect", predictedEffectColumn()),
                            col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                    .setHyperLink(hyperLink(new COSMICLinkExpression())).setStyle(linkStyle()),
                            col.column("Ploidy (TAF)", PatientDataSource.PLOIDY_TAF_FIELD)))
                        .setDataSource(PatientDataSource.fromVariants(report.variants(), reporterData)) :
                cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(
                cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                table,
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("*").setStyle(fontStyle()).setWidth(2),
                        cmp.text(geneMutationAddition).setStyle(fontStyle())),
                cmp.verticalGap(15),
                cmp.text("Implied Tumor Purity: " + report.impliedPurityString())
                        .setStyle(tableHeaderStyle()),
                cmp.verticalGap(15),
                cmp.text("Tumor Mutational Load: " + Integer.toString(report.mutationalLoad()) + " **")
                        .setStyle(tableHeaderStyle()),
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("**").setStyle(fontStyle()).setWidth(2),
                        cmp.text(mutationalLoadAddition).setStyle(fontStyle()))
        );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table;
        if (report.geneFusions().size() > 0) {
            table = cmp.subreport(
                        baseTable()
                            .fields(PatientDataSource.geneFusionFields())
                            .columns(
                                col.column("Gene A", PatientDataSource.GENE_FIELD).setFixedWidth(50),
                                col.column("Position A", PatientDataSource.POSITION_FIELD),
                                col.column("Gene A Context", PatientDataSource.SV_GENE_CONTEXT),
                                col.column("Gene B", PatientDataSource.SV_PARTNER_GENE_FIELD).setFixedWidth(50),
                                col.column("Position B", PatientDataSource.SV_PARTNER_POSITION_FIELD),
                                col.column("Gene B Context", PatientDataSource.SV_PARTNER_CONTEXT_FIELD),
                                col.column("SV Type", PatientDataSource.SV_TYPE_FIELD).setFixedWidth(30)
                            )
                            .setDataSource(PatientDataSource.fromGeneFusions(report.geneFusions()))
                    );
        } else {
            table = cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));
        }
        // @formatter:on

        return cmp.verticalList(cmp.text("Gene Fusions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table;
        if (report.geneDisruptions().size() > 0) {
            table = cmp.subreport(
                        baseTable()
                            .fields(PatientDataSource.geneDisruptionFields())
                            .columns(
                                col.column("Gene", PatientDataSource.GENE_FIELD).setFixedWidth(50),
                                col.column("Position", PatientDataSource.POSITION_FIELD),
                                col.column("Gene Context", PatientDataSource.SV_GENE_CONTEXT),
                                col.column("Orientation", PatientDataSource.SV_ORIENTATION_FIELD),
                                col.column("Partner", PatientDataSource.SV_PARTNER_POSITION_FIELD),
                                col.column("HGVS", PatientDataSource.SV_HGVS_FIELD),
                                col.column("Type", PatientDataSource.SV_TYPE_FIELD).setFixedWidth(30),
                                col.column("VAF", PatientDataSource.SV_VAF).setFixedWidth(30),
                                col.column("TAF", PatientDataSource.SV_TAF).setFixedWidth(30)
                            )
                            .setDataSource(PatientDataSource.fromGeneDisruptions(report.geneDisruptions()))
                    );
        } else {
            table = cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));
        }
        // @formatter:on

        return cmp.verticalList(cmp.text("Gene Disruptions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table = report.copyNumbers().size() > 0 ?
                cmp.subreport(baseTable().fields(PatientDataSource.copyNumberFields())
                        .columns(
                            col.column("Chromosome", PatientDataSource.CHROMOSOME_FIELD),
                            col.column("Band", PatientDataSource.BAND_FIELD),
                            col.column("Gene", PatientDataSource.GENE_FIELD),
                            col.column("Type", PatientDataSource.COPY_NUMBER_TYPE_FIELD),
                            col.column("Copies", PatientDataSource.COPY_NUMBER_FIELD))
                        .setDataSource(PatientDataSource.fromCopyNumbers(report.copyNumbers()))) :
                cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(
                cmp.text("Somatic Copy Numbers").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                table);
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> genePanelSection(@NotNull final HmfReporterData reporterData) {
        final long coverage = Math.round(reporterData.slicer().numberOfBases() / 1E6);
        final VerticalListBuilder section = toList("Details on the reported gene panel",
                Lists.newArrayList("The findings in this report are generated from whole-genome-sequencing analysis.",
                        "Findings are reported for the set of " + Integer.toString(reporterData.slicer().numberOfRegions())
                                + " genes (canonical transcripts) indicated below (covering " + coverage + " MBases)"));

        return section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP), createGenePanel(reporterData));
    }

    @NotNull
    private static ComponentBuilder<?, ?> createGenePanel(@NotNull final HmfReporterData reporterData) {
        //@formatter:off
        // KODU: Overwrite default font size to make the panel fit on one page.
        final int fontSize = 7;
        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle().setFontSize(fontSize)).fields(GenePanelDataSource.genePanelFields())
                    .columns(
                        col.emptyColumn().setFixedWidth(40),
                        col.column("Gene", GenePanelDataSource.GENE_FIELD).setFixedWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT_FIELD)
                                .setHyperLink(hyperLink(fieldTranscriptLink(GenePanelDataSource.TRANSCRIPT_FIELD)))
                                .setStyle(linkStyle().setFontSize(fontSize)).setFixedWidth(100),
                        col.column("Cosmic Type", GenePanelDataSource.TYPE_FIELD).setFixedWidth(75),
                        col.emptyColumn(),
                        col.column("Gene", GenePanelDataSource.GENE2_FIELD).setFixedWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT2_FIELD)
                                .setHyperLink(hyperLink(fieldTranscriptLink(GenePanelDataSource.TRANSCRIPT2_FIELD)))
                                .setStyle(linkStyle().setFontSize(fontSize)).setFixedWidth(100),
                        col.column("Cosmic Type", GenePanelDataSource.TYPE2_FIELD).setFixedWidth(75),
                        col.emptyColumn().setFixedWidth(40)))
                    .setDataSource(GenePanelDataSource.fromHmfReporterData(reporterData));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantFieldExplanationSection() {
        return toList("Details on reported genomic variant fields",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "The 'position' refers to the chromosome and start base of the variant with " + "respect to this reference genome.",
                        "The 'variant' displays what was expected as reference base and what " + "was found instead ('ref' > 'alt').",
                        "The 'depth (VAF)' displays the number of observations of the specific variant versus "
                                + "the total number of reads in this location in the format 'alt / total (%)'.",
                        "The 'predicted effect' provides additional information on the variant, including "
                                + "the change in coding sequence ('c.'), the change in protein ('p.') and "
                                + "the predicted impact on the final protein on the second line of this field.",
                        "The 'cosmic' fields display a link to the COSMIC database which contains "
                                + "additional information on the variant. If the variant could not be found in the "
                                + "COSMIC database, this field will be left blank. The Cosmic v76 database is used "
                                + "to look-up these IDs.",
                        "The implied tumor purity is the percentage of tumor DNA in the biopsy based on analysis of "
                                + "whole genome data.",
                        "The 'Ploidy (TAF)' field displays the tumor ploidy for the observed variant. The ploidy "
                                + "has been adjusted for the implied tumor purity (see above) and is shown as a "
                                + "proportion of A’s and B’s (e.g. AAABB for 3 copies A, and 2 copies B). "
                                + "The copy number is the sum of A’s and B’s. The TAF (Tumor adjusted Alternative "
                                + "Frequency) value refers to the alternative allele frequency after correction " + "for tumor purity.",
                        "The tumor mutational load is the total number of somatic missense variants found across"
                                + " the whole genome of the tumor biopsy."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberExplanationSection() {
        return toList("Details on reported copy numbers", Lists.newArrayList(
                "The lowest copy number value along the region of the canonical transcript is determined as "
                        + "a measure for the gene's copy number.",
                "Copy numbers are corrected for the implied tumor purity and represent the number of copies in the tumor DNA.",
                "Any gene with no copies is reported as loss.", "Any gene with at least 8 copies is reported as a gain.",
                "Any gene with more copies than 2.2 times the average tumor ploidy is reported as a gain."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> disclaimerSection() {
        return toList("Disclaimer", Lists.newArrayList("This test is not certified for diagnostic purposes.",
                "The findings in this report are not meant to be used for clinical decision making without validation of "
                        + "findings using certified assays.",
                "When no mutations are reported, the absence of mutations is not guaranteed."));
    }

    @NotNull
    private static VerticalListBuilder toList(@NotNull final String title, @NotNull final Iterable<String> lines) {
        final VerticalListBuilder list = cmp.verticalList();
        list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_HEADER_INDENT), cmp.text(title).setStyle(fontStyle().bold().setFontSize(11))),
                cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP));
        boolean isFirst = true;
        for (final String line : lines) {
            if (!isFirst) {
                list.add(cmp.verticalGap(DETAIL_TO_DETAIL_VERTICAL_GAP));
            }
            list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_DETAIL_INDENT), cmp.text("- ").setStyle(fontStyle()).setWidth(LIST_INDENT),
                    cmp.text(line).setStyle(fontStyle())));

            isFirst = false;
        }
        return list;
    }

    @NotNull
    private static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
    }

    @NotNull
    private static StyleBuilder sectionHeaderStyle() {
        return fontStyle().bold().setFontSize(12).setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);
    }

    @NotNull
    private static StyleBuilder tableHeaderStyle() {
        return fontStyle().bold()
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setFontSize(10)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataTableStyle() {
        return dataStyle().setBorder(stl.pen1Point());
    }

    @NotNull
    private static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }

    @NotNull
    private static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
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

    @NotNull
    private static ComponentBuilder<?, ?> predictedEffectColumn() {
        return cmp.verticalList(cmp.horizontalList(cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.CONSEQUENCE_FIELD))).setFixedWidth(170);
    }
}
