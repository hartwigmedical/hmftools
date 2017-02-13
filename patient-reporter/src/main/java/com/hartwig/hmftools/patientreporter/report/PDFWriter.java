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
import java.util.Collection;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter {

    private static final String FONT = "Times New Roman";

    private static final int TEXT_HEADER_INDENT = 30;
    private static final int TEXT_DETAIL_INDENT = 40;
    private static final int LIST_INDENT = 5;
    private static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;
    private static final int DETAIL_TO_DETAIL_VERTICAL_GAP = 4;
    private static final int SECTION_VERTICAL_GAP = 25;

    @NotNull
    private final String outputDirectory;
    @NotNull
    private final String hmfLogo;
    @NotNull
    private final Slicer hmfSlicingRegion;

    public PDFWriter(@NotNull final String outputDirectory, @NotNull final String hmfLogo,
            @NotNull final Slicer hmfSlicingRegion) {
        this.outputDirectory = outputDirectory;
        this.hmfLogo = hmfLogo;
        this.hmfSlicingRegion = hmfSlicingRegion;
    }

    @NotNull
    public String writeSequenceReport(@NotNull final PatientReport report) throws FileNotFoundException, DRException {
        final String fileName = fileName(report.sample());
        final JasperReportBuilder jasperReportBuilder = generatePatientReport(report, hmfLogo, hmfSlicingRegion);

        jasperReportBuilder.toPdf(new FileOutputStream(fileName));

        return fileName;
    }

    @NotNull
    public String writeNonSequenceableReport(@NotNull final String sample, @NotNull String tumorType,
            @NotNull final NotSequenceableReason reason) throws FileNotFoundException, DRException {
        final String fileName = fileName(sample);
        final JasperReportBuilder jasperReportBuilder = generateNotSequenceableReport(sample, tumorType, reason,
                hmfLogo);

        jasperReportBuilder.toPdf(new FileOutputStream(fileName));

        return fileName;
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return outputDirectory + File.separator + sample + "_hmf_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generateNotSequenceableReport(@NotNull final String sample,
            @NotNull final String tumorType, @NotNull final NotSequenceableReason reason,
            @NotNull final String hmfLogoPath) {
        // @formatter:off
        final ComponentBuilder<?, ?> report =
                cmp.verticalList(
                        mainPageTopSection(sample, tumorType, hmfLogoPath),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageNotSequenceableSection(reason));
        // @formatter:on

        return report().noData(report);
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final String hmfLogoPath, @NotNull final Slicer hmfSlicingRegion) {
        // @formatter:off
        final ComponentBuilder<?, ?> reportMainPage =
                cmp.verticalList(
                        mainPageTopSection(report.sample(), report.tumorType(), hmfLogoPath),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        mainPageAboutSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantReport(report),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberReport(report));

        final ComponentBuilder<?, ?> helpPage =
                cmp.verticalList(
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        cmp.text("HMF Sequencing Report - Additional Information").setStyle(sectionHeaderStyle()),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        filteringSection(hmfSlicingRegion),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        variantFieldExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        copyNumberExplanationSection(),
                        cmp.verticalGap(SECTION_VERTICAL_GAP),
                        disclaimerSection()
                );

        final ComponentBuilder<?, ?> totalReport =
                cmp.multiPageList().add(reportMainPage).newPage().add(helpPage);
        // @formatter:on

        return report().noData(totalReport);
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageTopSection(@NotNull final String sample,
            @NotNull final String tumorType, @NotNull final String hmfLogoPath) {
        // @formatter:off
        final ComponentBuilder<?, ?> mainDiagnosisInfo = cmp.horizontalList(
                cmp.verticalList(
                        cmp.text("Report Date").setStyle(tableHeaderStyle()),
                        cmp.currentDate().setPattern("dd-MMM-yyyy").setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Tumor Type").setStyle(tableHeaderStyle()),
                        cmp.text(tumorType).setStyle(dataTableStyle()))
        );

        return cmp.horizontalList(
                cmp.image(hmfLogoPath),
                cmp.verticalList(
                        cmp.text("HMF Sequencing Report - " + sample).
                                setStyle(fontStyle().bold().setFontSize(14)
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
                "For additional questions, please contact us via info@hartwigmedicalfoundation.nl."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> mainPageNotSequenceableSection(@NotNull final NotSequenceableReason reason) {

        // @formatter:off
        return cmp.verticalList(
                cmp.text("Notification of inadequate tumor sample:").setStyle(tableHeaderStyle().setFontSize(12))
                        .setHeight(20),
                cmp.text("Insufficient percentage of tumor cells").setStyle(dataTableStyle().setFontSize(12))
                        .setHeight(20),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The received tumor sample for this patient was inadequate to obtain a reliable sequencing " +
                                "result. For sequencing we require a minimum of 30% tumor cells. " +
                                "Therefore whole genome sequencing cannot be performed, unless additional fresh tumor " +
                                "material can be provided for a new assessment.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("When possible, please resubmit using the same CPCT-number. In case additional tumor " +
                        "material cannot be provided, please be notified that the patient will not be evaluable " +
                        "for the CPCT-02 study.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("For questions, please contact us via info@hartwigmedicalfoundation.nl")
                        .setStyle(fontStyle()));

        // @formatter:on
        //
        //        final String reasonString;
        //        switch (reason) {
        //            case LOW_DNA_YIELD:
        //                reasonString = "The biopsy has not been analyzed because of low DNA yield from the biopsy";
        //                break;
        //            case LOW_TUMOR_PERCENTAGE:
        //                reasonString = "The biopsy could not be analyzed because the biopsy contained less than 30% tumor cells";
        //                break;
        //            default:
        //                reasonString = "ERROR";
        //        }
        //
        //        return toList("About this report", Lists.newArrayList(reasonString,
        //                "For additional questions, please contact us via info@hartwigmedicalfoundation.nl."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table = report.variants().size() > 0 ?
                cmp.subreport(baseTable().fields(PatientDataSource.variantFields())
                        .columns(
                            col.column("Gene", PatientDataSource.GENE_FIELD),
                            col.column("Position", PatientDataSource.POSITION_FIELD),
                            col.column("Variant", PatientDataSource.VARIANT_FIELD),
                            transcriptColumn(),
                            col.componentColumn("Predicted Effect", predictedEffectColumn()),
                            col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                    .setHyperLink(hyperLink(new COSMICLinkExpression())).setStyle(linkStyle()),
                            col.column("VAF", PatientDataSource.ALLELE_FREQUENCY_FIELD)))
                        .setDataSource(PatientDataSource.fromVariants(report.variants())) :
                cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(
                cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                table,
                cmp.verticalGap(15),
                cmp.text("Tumor Mutational Load: " + Integer.toString(report.mutationalLoad()))
                        .setStyle(tableHeaderStyle())
        );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberReport(@NotNull final PatientReport report) {
        // @formatter:off
        final ComponentBuilder<?, ?> table = report.copyNumbers().size() > 0 ?
                cmp.subreport(baseTable().fields(PatientDataSource.copyNumberFields())
                        .columns(
                            col.column("Gene", PatientDataSource.GENE_FIELD),
                            transcriptColumn(),
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
    private static ComponentBuilder<?, ?> filteringSection(@NotNull final Slicer hmfSlicingRegion) {
        final long coverage = Math.round(hmfSlicingRegion.numberOfBases() / 1E6);
        final VerticalListBuilder section = toList("Details on filtering",
                Lists.newArrayList("The findings in this report are generated from whole-genome-sequencing analysis.",
                        "The results are filtered on the canonical transcripts for the set of below "
                                + Integer.toString(hmfSlicingRegion.numberOfRegions()) + " genes (covering " + coverage
                                + " MBases)",
                        "The definition of canonical transcripts can be found on http://www.ensembl.org/Help/Glossary?id=346"));

        return section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                createGenePanel(hmfSlicingRegion.regions()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> createGenePanel(@NotNull final Collection<GenomeRegion> regions) {
        final Collection<String> genes = Sets.newTreeSet();
        for (final GenomeRegion region : regions) {
            final HMFSlicingAnnotation annotation = HMFSlicingAnnotation.fromGenomeRegion(region);
            // KODU: The annotation should always be present on the HMF slicing regions!
            assert annotation != null;
            genes.add(annotation.gene());
        }
        final VerticalListBuilder table = cmp.verticalList();
        final int nrOfGenesPerRow = 10;

        long nrOfRowsNeeded = Math.round((double) genes.size() / nrOfGenesPerRow);
        nrOfRowsNeeded = (nrOfRowsNeeded * nrOfGenesPerRow < genes.size()) ? nrOfRowsNeeded + 1 : nrOfRowsNeeded;

        for (int i = 0; i < nrOfRowsNeeded; i++) {
            final HorizontalListBuilder row = cmp.horizontalList();
            for (int j = 0; j < nrOfGenesPerRow; j++) {
                int index = i * nrOfGenesPerRow + j + 1;
                final String gene = index > genes.size() ? Strings.EMPTY : (String) genes.toArray()[index - 1];
                row.add(cmp.text(gene).setStyle(dataTableStyle()));
            }
            table.add(row);
        }
        return table;
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantFieldExplanationSection() {
        return toList("Details on reported genomic variant fields",
                Lists.newArrayList("The analysis is based on reference genome version GRCh37.",
                        "The 'position' refers to the chromosome and start base of the variant with "
                                + "respect to this reference genome.",
                        "The 'variant' displays what was expected as reference base and what "
                                + "was found instead ('ref' > 'alt').",
                        "The 'transcript' provides a link to the ensembl definition of the transcript "
                                + "used for filtering.",
                        "The 'predicted effect' provides additional information on the variant, including "
                                + "the change in coding sequence ('c.'), the change in amino acid ('a.') and "
                                + "the predicted impact on the final protein on the second line of this field.",
                        "The 'cosmic' fields display a link to the COSMIC database which contains "
                                + "additional information on the variant. If the variant could not be found in the "
                                + "COSMIC database, this field will be left blank. The Cosmic v76 database is used "
                                + "to look-up these IDs.",
                        "The 'VAF' fields displays the variant allele frequency. The first number is "
                                + "the number of observations of the variant, and the second number is the total "
                                + "number of observations on this position. The number within parentheses is the "
                                + "allele frequency (the two numbers divided by each other).",
                        "The tumor mutational load is the total number of somatic missense variants found across "
                                + " the whole genome of the tumor biopsy."));
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberExplanationSection() {
        return toList("Details on reported copy numbers",
                Lists.newArrayList("Copy numbers are determined for the genes filtered for in this report.",
                        "The lowest copy number value along the region of the canonical transcript is determined as a measure for the gene's copy number.",
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
                new Color(210, 210, 210));
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8).setHorizontalTextAlignment(
                HorizontalTextAlignment.CENTER).setVerticalTextAlignment(VerticalTextAlignment.MIDDLE);
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
    private static TextColumnBuilder<?> transcriptColumn() {
        return col.column("Transcript", PatientDataSource.TRANSCRIPT_FIELD).setWidth(150).setHyperLink(
                hyperLink(new TranscriptLinkExpression())).setStyle(linkStyle());
    }

    @NotNull
    private static ComponentBuilder<?, ?> predictedEffectColumn() {
        return cmp.verticalList(
                cmp.horizontalList(cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.EFFECT_FIELD)));
    }
}
