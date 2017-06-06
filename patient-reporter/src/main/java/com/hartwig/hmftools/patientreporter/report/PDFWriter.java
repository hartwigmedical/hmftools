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
import java.nio.file.Files;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.genePanel.GenePanelModel;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotationFactory;

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
    private static final String VERSION = "2.0";

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
    private final String reportLogo;

    public PDFWriter(@NotNull final String reportDirectory, @NotNull final String reportLogo) {
        this.reportDirectory = reportDirectory;
        this.reportLogo = reportLogo;
    }

    @NotNull
    @Override
    public String writeSequenceReport(@NotNull final PatientReport report, @NotNull final Slicer hmfSlicingRegion,
            @NotNull final DrupFilter drupFilter, @NotNull final GenePanelModel genePanelModel)
            throws FileNotFoundException, DRException {
        final JasperReportBuilder reportBuilder = generatePatientReport(report, reportLogo, hmfSlicingRegion,
                drupFilter);

        return writeReport(report.sample(), reportBuilder);
    }

    @NotNull
    @Override
    public String writeNonSequenceableReport(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage, @NotNull final NotSequenceableReason reason)
            throws FileNotFoundException, DRException {
        final JasperReportBuilder reportBuilder = generateNotSequenceableReport(sample, tumorType, tumorPercentage,
                reason, reportLogo);

        return writeReport(sample, reportBuilder);
    }

    @NotNull
    private String writeReport(@NotNull final String sample, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        final String fileName = fileName(sample);
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write report as it already exists: " + fileName);
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created patient report at " + fileName);
        }
        return fileName;
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_hmf_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generateNotSequenceableReport(@NotNull final String sample,
            @NotNull final String tumorType, @NotNull final String tumorPercentageString,
            @NotNull final NotSequenceableReason reason, @NotNull final String reportLogoPath) {
        // @formatter:off
        final ComponentBuilder<?, ?> report =
                cmp.verticalList(
                        mainPageTopSection(sample, tumorType, tumorPercentageString, reportLogoPath),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageNotSequenceableSection(reason));
        // @formatter:on

        return report().noData(report);
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final String reportLogoPath, @NotNull final Slicer hmfSlicingRegion,
            @NotNull final DrupFilter drupFilter) {
        // @formatter:off
        final ComponentBuilder<?, ?> reportMainPage =
                cmp.verticalList(
                        mainPageTopSection(report.sample(), report.tumorType(), report.tumorPercentageString(),
                                reportLogoPath),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageAboutSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantReport(report, drupFilter),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberReport(report));

        final ComponentBuilder<?, ?> genePanelPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report v" + VERSION + " - Gene Panel Information")
                                .setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        genePanelSection(hmfSlicingRegion)
                );

        final ComponentBuilder<?, ?> additionalInfoPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report v" + VERSION + " - Additional Information")
                                .setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantFieldExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        disclaimerSection()
                );

        final ComponentBuilder<?, ?> totalReport =
                cmp.multiPageList().add(reportMainPage).newPage().add(genePanelPage).newPage().add(additionalInfoPage);
        // @formatter:on

        return report().noData(totalReport);
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageTopSection(@NotNull final String sample,
            @NotNull final String tumorType, @NotNull final String tumorPercentage,
            @NotNull final String reportLogoPath) {
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
                cmp.image(reportLogoPath),
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
        return toList("About this report", Lists.newArrayList(
                "This test is performed for research purpose and is not meant to be used for "
                        + "clinical decision making without further validation of findings.",
                "Additional information on the various fields can be found on the final page of this report.",
                "For DRUP-specific questions, please contact the DRUP study team at DRUP@nki.nl.",
                "For other questions, please contact us via info@hartwigmedicalfoundation.nl."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageNotSequenceableSection(@NotNull final NotSequenceableReason reason) {
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
                cmp.text("When possible, please resubmit using the same CPCT-number. In case additional tumor " +
                        "material cannot be provided, please be notified that the patient will not be evaluable " +
                        "for the CPCT-02 study.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("For questions, please contact us via info@hartwigmedicalfoundation.nl")
                        .setStyle(fontStyle()));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantReport(@NotNull final PatientReport report,
            @NotNull final DrupFilter drupFilter) {
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
                            col.column("Depth", PatientDataSource.READ_DEPTH_FIELD),
                            col.componentColumn("Predicted Effect", predictedEffectColumn()),
                            col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                    .setHyperLink(hyperLink(new COSMICLinkExpression())).setStyle(linkStyle()),
                            col.column("BAF (VAF)", PatientDataSource.BAF_VAF_FIELD)))
                        .setDataSource(PatientDataSource.fromVariants(report.variants(), drupFilter)) :
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
    private static ComponentBuilder<?, ?> copyNumberReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table = report.copyNumbers().size() > 0 ?
                cmp.subreport(baseTable().fields(PatientDataSource.copyNumberFields())
                        .columns(
                            col.column("Chromosome", PatientDataSource.CHROMOSOME_FIELD),
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
    private static ComponentBuilder<?, ?> genePanelSection(@NotNull final Slicer hmfSlicingRegion) {
        final long coverage = Math.round(hmfSlicingRegion.numberOfBases() / 1E6);
        final VerticalListBuilder section = toList("Details on gene panel",
                Lists.newArrayList("The findings in this report are generated from whole-genome-sequencing analysis.",
                        "The results are filtered on the below canonical transcripts for a set of "
                                + Integer.toString(hmfSlicingRegion.numberOfRegions()) + " genes (covering " + coverage
                                + " MBases)",
                        "The definition of canonical transcripts can be found on http://www.ensembl.org/Help/Glossary?id=346"));

        return section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                createGenePanel(hmfSlicingRegion.regions()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> createGenePanel(@NotNull final Collection<GenomeRegion> regions) {
        //@formatter:off
        final List<HMFSlicingAnnotation> annotations = Lists.newArrayList();
        for (final GenomeRegion region : regions) {
            final HMFSlicingAnnotation annotation = HMFSlicingAnnotationFactory.fromGenomeRegion(region);
            // KODU: The annotation should always be present on the HMF slicing regions!
            assert annotation != null;
            annotations.add(annotation);
        }
        annotations.sort(Comparator.comparing(HMFSlicingAnnotation::gene));

        return cmp.subreport(
                baseTable().fields(GenePanelDataSource.genePanelFields())
                    .columns(
                        col.emptyColumn().setFixedWidth(40),
                        col.column("Gene", GenePanelDataSource.GENE_FIELD).setFixedWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT_FIELD)
                                .setHyperLink(hyperLink(fieldTranscriptLink(GenePanelDataSource.TRANSCRIPT_FIELD)))
                                .setStyle(linkStyle()).setFixedWidth(100),
                        col.emptyColumn(),
                        col.column("Gene", GenePanelDataSource.GENE2_FIELD).setFixedWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT2_FIELD)
                                .setHyperLink(hyperLink(fieldTranscriptLink(GenePanelDataSource.TRANSCRIPT2_FIELD)))
                                .setStyle(linkStyle()).setFixedWidth(100),
                        col.emptyColumn(),
                        col.column("Gene", GenePanelDataSource.GENE3_FIELD).setFixedWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT3_FIELD)
                                .setHyperLink(hyperLink(fieldTranscriptLink(GenePanelDataSource.TRANSCRIPT3_FIELD)))
                                .setStyle(linkStyle()).setFixedWidth(100),
                        col.emptyColumn().setFixedWidth(40)))
                    .setDataSource(GenePanelDataSource.fromHMFSlicingAnnotations(annotations));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantFieldExplanationSection() {
        return toList("Details on reported genomic variant fields",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "The 'position' refers to the chromosome and start base of the variant with "
                                + "respect to this reference genome.",
                        "The 'variant' displays what was expected as reference base and what "
                                + "was found instead ('ref' > 'alt').",
                        "The 'read depth' displays the number of observations of the specific variant versus "
                                + "the number of reference reads in this location in the format 'alt / ref (%)'.",
                        "The 'predicted effect' provides additional information on the variant, including "
                                + "the change in coding sequence ('c.'), the change in protein ('p.') and "
                                + "the predicted impact on the final protein on the second line of this field.",
                        "The 'cosmic' fields display a link to the COSMIC database which contains "
                                + "additional information on the variant. If the variant could not be found in the "
                                + "COSMIC database, this field will be left blank. The Cosmic v76 database is used "
                                + "to look-up these IDs.",
                        "The 'BAF (VAF)' fields displays the purity-adjusted beta allele frequency of this location "
                                + "as a proportion of A's and B's. The copy number is the sum of A's and B's. The "
                                + "VAF in parentheses displays the allele frequency of this variant in the tumor "
                                + "corrected for purity.",
                        "The implied tumor purity is the percentage of tumor cells in the biopsy based on analysis of "
                                + " whole genome data.",
                        "The tumor mutational load is the total number of somatic missense variants found across "
                                + " the whole genome of the tumor biopsy."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberExplanationSection() {
        return toList("Details on reported copy numbers",
                Lists.newArrayList("Copy numbers are determined for the genes filtered for in this report.",
                        "The lowest copy number value along the region of the canonical transcript is determined as "
                                + "a measure for the gene's copy number.",
                        "Any gene with a value of 0 (loss) or >3 (gain) is included in the list of findings."));
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
        list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_HEADER_INDENT),
                cmp.text(title).setStyle(fontStyle().bold().setFontSize(11))),
                cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP));
        boolean isFirst = true;
        for (final String line : lines) {
            if (!isFirst) {
                list.add(cmp.verticalGap(DETAIL_TO_DETAIL_VERTICAL_GAP));
            }
            list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_DETAIL_INDENT),
                    cmp.text("- ").setStyle(fontStyle()).setWidth(LIST_INDENT), cmp.text(line).setStyle(fontStyle())));

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
        return fontStyle().bold().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER).setVerticalTextAlignment(
                VerticalTextAlignment.MIDDLE).setFontSize(10).setBorder(stl.pen1Point()).setBackgroundColor(
                BORKIE_COLOR).setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8).setHorizontalTextAlignment(
                HorizontalTextAlignment.CENTER).setVerticalTextAlignment(VerticalTextAlignment.MIDDLE).setPadding(
                PADDING);
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
                return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(
                        field.getName());
            }
        };
    }

    @NotNull
    private static ComponentBuilder<?, ?> predictedEffectColumn() {
        return cmp.verticalList(
                cmp.horizontalList(cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.EFFECT_FIELD))).setFixedWidth(170);
    }
}
