package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.monospaceBaseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.report.components.GenePanelSection;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.VariantDataSource;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableCircosPage;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableExplanationPage;
import com.hartwig.hmftools.patientreporter.report.pages.NonSequenceablePage;
import com.hartwig.hmftools.patientreporter.report.pages.SampleDetailsPage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(PDFWriter.class);

    @NotNull
    private final String reportDirectory;

    public PDFWriter(@NotNull final String reportDirectory) {
        this.reportDirectory = reportDirectory;
    }

    @Override
    public void writeSequenceReport(@NotNull final SequencedPatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generatePatientReport(report, reporterData);
        writeReport(fileName(report.sampleReport().sampleId(), "_hmf_report.pdf"), reportBuilder);
        final JasperReportBuilder evidenceReportBuilder = EvidenceReport.generate(report);
        writeReport(fileName(report.sampleReport().sampleId(), "_evidence_items.pdf"), evidenceReportBuilder);
    }

    @Override
    public void writeNonSequenceableReport(@NotNull final NotSequencedPatientReport report) throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generateNotSequenceableReport(report);
        writeReport(fileName(report.sampleReport().sampleId(), "_hmf_report.pdf"), reportBuilder);
    }

    private void writeReport(@NotNull final String fileName, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write " + fileName + " as it already exists.");
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created patient report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample, @NotNull final String suffix) {
        return reportDirectory + File.separator + sample + suffix;
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generateNotSequenceableReport(@NotNull final NotSequencedPatientReport report) throws IOException {
        // MIVO: hack to get page footers working; the footer band and noData bands are exclusive, see additional comment below for details
        final DRDataSource singleItemDataSource = new DRDataSource("item");
        singleItemDataSource.add(new Object());

        return report().pageFooter(cmp.pageXslashY())
                .lastPageFooter(cmp.verticalList(signatureFooter(report.signaturePath()),
                        cmp.pageXslashY(),
                        cmp.text("End of report.").setStyle(stl.style().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .addDetail(NonSequenceablePage.of(report).reportComponent())
                .setDataSource(singleItemDataSource);
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final SequencedPatientReport report,
            @NotNull final HmfReporterData reporterData) throws IOException {
        final ComponentBuilder<?, ?> reportMainPage = cmp.verticalList(MainPageTopSection.buildWithImpliedPurity("HMF Sequencing Report",
                report.sampleReport(),
                report.impliedPurityString()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                mainPageAboutSection(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                snvIndelReport(report, reporterData),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                copyNumberReport(report),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneDisruptionReport(report),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneFusionReport(report));

        final ComponentBuilder<?, ?> genePanelPage = cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(Commons.TITLE + " - Gene Panel Information").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                GenePanelSection.build(reporterData));

        final ComponentBuilder<?, ?> totalReport = cmp.multiPageList()
                .add(reportMainPage)
                .newPage()
                .add(ImmutableCircosPage.of(report.circosPath()).reportComponent())
                .newPage()
                .add(genePanelPage)
                .newPage()
                .add(ImmutableExplanationPage.builder().build().reportComponent())
                .newPage()
                .add(SampleDetailsPage.of(report).reportComponent());

        // MIVO: hack to get page footers working; the footer band and noData bands are exclusive:
        //  - footerBand, detailBand, etc are shown when data source is not empty
        //  - noData band is shown when data source is empty; intended to be used when there is no data to show in the report
        //  (e.g. would be appropriate to be used for notSequenceableReport)
        //
        // more info: http://www.dynamicreports.org/examples/bandreport
        //
        // todo: fix when implementing new report layout

        final DRDataSource singleItemDataSource = new DRDataSource("item");
        singleItemDataSource.add(new Object());

        return report().pageFooter(cmp.pageXslashY())
                .lastPageFooter(cmp.verticalList(signatureFooter(report.signaturePath()),
                        cmp.pageXslashY(),
                        cmp.text("End of report.").setStyle(stl.style().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .addDetail(totalReport)
                .setDataSource(singleItemDataSource);
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageAboutSection() {
        return toList("About this report",
                Lists.newArrayList("This report is based on tests that are performed under ISO/ICE-17025:2005 accreditation.",
                        "For DRUP-specific questions, please contact the DRUP study team at DRUP@nki.nl.",
                        "For other questions, please contact us via info@hartwigmedicalfoundation.nl."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> snvIndelReport(@NotNull final SequencedPatientReport report,
            @NotNull final HmfReporterData reporterData) {
        final String mutationalLoadAddition =
                "Patients with a mutational load over 140 could be " + "eligible for immunotherapy within DRUP.";

        final String geneMutationAddition = "Marked genes (*) are included in the DRUP study and indicate potential "
                + "eligibility in DRUP. Please note that the marking is NOT based on the specific variant reported for "
                + "this sample, but only on a gene-level.";

        final ComponentBuilder<?, ?> table =
                report.variants().size() > 0
                        ? cmp.subreport(monospaceBaseTable().fields(VariantDataSource.variantFields())
                        .columns(col.column("Gene", VariantDataSource.GENE_FIELD).setFixedWidth(50),
                                col.column("Position", VariantDataSource.POSITION_FIELD),
                                col.column("Variant", VariantDataSource.VARIANT_FIELD),
                                col.column("Depth (VAF)", VariantDataSource.DEPTH_VAF_FIELD),
                                col.componentColumn("Predicted Effect", predictedEffectColumn()),
                                col.column("Cosmic", VariantDataSource.COSMIC_FIELD)
                                        .setHyperLink(hyperLink(new COSMICLinkExpression()))
                                        .setStyle(linkStyle()),
                                col.column("Ploidy (TAF)", VariantDataSource.PLOIDY_TAF_FIELD)))
                        .setDataSource(VariantDataSource.fromVariants(report.variants(), reporterData))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                table,
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("*").setStyle(fontStyle()).setWidth(2),
                        cmp.text(geneMutationAddition).setStyle(fontStyle())),
                cmp.verticalGap(15),
                cmp.text("Tumor Mutational Load: " + Integer.toString(report.mutationalLoad()) + " **").setStyle(tableHeaderStyle()),
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("**").setStyle(fontStyle()).setWidth(2),
                        cmp.text(mutationalLoadAddition).setStyle(fontStyle())));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                report.geneCopyNumbers().size() > 0
                        ? cmp.subreport(monospaceBaseTable().fields(GeneCopyNumberDataSource.copyNumberFields())
                        .columns(col.column("Position", GeneCopyNumberDataSource.POSITION_FIELD),
                                col.column("Gene", GeneCopyNumberDataSource.GENE_FIELD),
                                col.column("Type", GeneCopyNumberDataSource.GAIN_OR_LOSS_FIELD),
                                col.column("Copies", GeneCopyNumberDataSource.COPY_NUMBER_FIELD))
                        .setDataSource(GeneCopyNumberDataSource.fromCopyNumbers(report.geneCopyNumbers())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Copy Numbers").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.geneDisruptions().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(GeneDisruptionDataSource.geneDisruptionFields())
                .columns(col.column("Gene", GeneDisruptionDataSource.GENE_FIELD).setFixedWidth(50),
                        col.column("Position", GeneDisruptionDataSource.POSITION_FIELD),
                        col.column("Gene Context", GeneDisruptionDataSource.SV_GENE_CONTEXT),
                        col.column("Orientation", GeneDisruptionDataSource.SV_ORIENTATION_FIELD),
                        col.column("Partner", GeneDisruptionDataSource.SV_PARTNER_POSITION_FIELD),
                        col.column("Type", GeneDisruptionDataSource.SV_TYPE_FIELD).setFixedWidth(30),
                        col.column("VAF", GeneDisruptionDataSource.SV_VAF).setFixedWidth(30))
                .setDataSource(GeneDisruptionDataSource.fromGeneDisruptions(report.geneDisruptions())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Disruptions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                report.geneFusions().size() > 0
                        ? cmp.subreport(monospaceBaseTable().fields(GeneFusionDataSource.geneFusionFields())
                        .columns(col.column("5' Gene", GeneFusionDataSource.GENE_FIELD).setFixedWidth(50),
                                col.column("5' Position", GeneFusionDataSource.POSITION_FIELD),
                                col.column("5' Gene Context", GeneFusionDataSource.SV_GENE_CONTEXT),
                                col.column("3' Gene", GeneFusionDataSource.SV_PARTNER_GENE_FIELD).setFixedWidth(50),
                                col.column("3' Position", GeneFusionDataSource.SV_PARTNER_POSITION_FIELD),
                                col.column("3' Gene Context", GeneFusionDataSource.SV_PARTNER_CONTEXT_FIELD),
                                col.column("SV Type", GeneFusionDataSource.SV_TYPE_FIELD).setFixedWidth(30),
                                col.column("VAF", GeneFusionDataSource.SV_VAF).setFixedWidth(30))
                        .setDataSource(GeneFusionDataSource.fromGeneFusions(report.geneFusions())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Fusions").setStyle(sectionHeaderStyle()), cmp.verticalGap(6), table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> signatureFooter(@NotNull final String signaturePath) {
        // @formatter:off
        return cmp.horizontalList(
                cmp.horizontalGap(370),
                cmp.xyList()
                    .add(40, 5, cmp.image(signaturePath))
                    .add(0, 0,cmp.text("Edwin Cuppen,"))
                    .add(0, 15, cmp.text("Director Hartwig Medical Foundation").setWidth(190)),
                cmp.horizontalGap(10));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> predictedEffectColumn() {
        return cmp.verticalList(cmp.horizontalList(cmp.text(DataExpression.fromField(VariantDataSource.HGVS_CODING_FIELD)),
                cmp.text(DataExpression.fromField(VariantDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(VariantDataSource.CONSEQUENCE_FIELD))).setFixedWidth(170);
    }
}
