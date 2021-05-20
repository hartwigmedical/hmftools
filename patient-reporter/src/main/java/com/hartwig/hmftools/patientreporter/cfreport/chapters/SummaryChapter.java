package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.hartwig.hmftools.patientreporter.cfreport.data.GainsAndLosses;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.HomozygousDisruptions;
import com.hartwig.hmftools.patientreporter.cfreport.data.Pharmacogenetics;
import com.hartwig.hmftools.patientreporter.cfreport.data.SomaticVariants;
import com.hartwig.hmftools.patientreporter.cfreport.data.TumorPurity;
import com.hartwig.hmftools.patientreporter.cfreport.data.VirusBreakends;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SummaryChapter implements ReportChapter {

    private static final float TABLE_SPACER_HEIGHT = 5;
    private static final DecimalFormat SINGLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat DOUBLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.##");

    @NotNull
    private final AnalysedPatientReport patientReport;

    public SummaryChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String name() {
        if (patientReport.isCorrectedReport()) {
            return "DNA Analysis Report (Corrected)";
        } else {
            if (patientReport.qsFormNumber().equals(QsFormNumber.FOR_209.display())) {
                return "DNA Analysis Report - Low Sensitivity";
            } else {
                return "DNA Analysis Report";
            }
        }
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public boolean hasCompleteSidebar() {
        return true;
    }

    private GenomicAnalysis analysis() {
        return patientReport.genomicAnalysis();
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        reportDocument.add(new Paragraph("Summary").addStyle(ReportResources.chapterTitleStyle()));
        reportDocument.add(TumorLocationAndTypeTable.createBiopsyLocationAndTumorLocation(patientReport.sampleReport()
                .primaryTumorLocationString(), patientReport.sampleReport().biopsyLocationString(), contentWidth()));
        reportDocument.add(new Paragraph());
        reportDocument.add(TumorLocationAndTypeTable.createTumorType(patientReport.sampleReport().primaryTumorTypeString(),
                contentWidth()));
        reportDocument.add(new Paragraph("\nThe information regarding 'primary tumor location', 'primary tumor type' and 'biopsy location'"
                + "  \nis based on information received from the originating hospital.").addStyle(ReportResources.subTextStyle()));

        renderClinicalConclusionText(reportDocument);
        renderTreatmentIndications(reportDocument);
        renderTumorCharacteristics(reportDocument);
        renderGenomicAlterations(reportDocument);
        renderPeach(reportDocument);
    }

    private void renderClinicalConclusionText(@NotNull Document reportDocument) {
        String text = patientReport.clinicalSummary();
        if (text.isEmpty()) {
            if (!analysis().hasReliablePurity()) {
                text = "Of note, WGS analysis indicated a very low abundance of genomic aberrations, which can be caused "
                        + "by a low tumor percentage in the received tumor material or due to genomic very stable/normal tumor type. "
                        + "As a consequence no reliable tumor purity assessment is possible and no information regarding "
                        + "mutation copy number and tVAF can be provided.";
            } else if (analysis().impliedPurity() < ReportResources.PURITY_CUTOFF) {
                double impliedPurityPercentage =
                        MathUtil.mapPercentage(analysis().impliedPurity(), TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
                text = "Due to the lower sensitivity (" + DataUtil.formatPercentage(impliedPurityPercentage) + ") "
                        + "of this test potential (subclonal) DNA aberrations might not have been detected using this test. " + ""
                        + "This result should therefore be considered with caution.";
            }
        }

        if (!text.isEmpty()) {
            Div div = createSectionStartDiv(contentWidth());
            div.add(new Paragraph("Clinical Conclusion").addStyle(ReportResources.sectionTitleStyle()));

            div.add(new Paragraph(text).setWidth(contentWidth()).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(11));

            reportDocument.add(div);
        }
    }

    private void renderTreatmentIndications(@NotNull Document reportDocument) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCellSummary()
                .add(new Paragraph("Treatment options (tumor-type specific)").addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableUtil.createLayoutCell(4, 2).setHeight(TABLE_SPACER_HEIGHT));

        int therapyEventCount = EvidenceItems.uniqueEventCount(analysis().tumorSpecificEvidence());
        int therapyCount = EvidenceItems.uniqueTherapyCount(analysis().tumorSpecificEvidence());
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with therapy indication").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(therapyEventCount, therapyCount, "treatment(s)"));

        int trialEventCount = ClinicalTrials.uniqueEventCount(analysis().clinicalTrials());
        int trialCount = ClinicalTrials.uniqueTrialCount(analysis().clinicalTrials());
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with clinical trial eligibility").addStyle(
                ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(trialEventCount, trialCount, "trial(s)"));

        div.add(table);

        reportDocument.add(div);
    }

    private void renderTumorCharacteristics(@NotNull Document reportDocument) {
        boolean hasReliablePurity = analysis().hasReliablePurity();

        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, .33f, .66f }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Tumor characteristics").setVerticalAlignment(VerticalAlignment.TOP)
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT));

        double impliedPurity = analysis().impliedPurity();
        double impliedPurityPercentage = MathUtil.mapPercentage(impliedPurity, TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
        renderTumorPurity(hasReliablePurity,
                DataUtil.formatPercentage(impliedPurityPercentage),
                impliedPurity,
                TumorPurity.RANGE_MIN,
                TumorPurity.RANGE_MAX,
                table);

        String molecularTissuePrediction = hasReliablePurity && patientReport.qsFormNumber().equals(QsFormNumber.FOR_080.display())
                ? patientReport.molecularTissueOrigin().conclusion()
                : DataUtil.NA_STRING;
        Style dataStyleMolecularTissuePrediction = hasReliablePurity && patientReport.qsFormNumber().equals(QsFormNumber.FOR_080.display())
                ? ReportResources.dataHighlightStyle()
                : ReportResources.dataHighlightNaStyle();

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Molecular tissue of origin prediction").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(molecularTissuePrediction).addStyle(
                dataStyleMolecularTissuePrediction)));

        Style dataStyle = hasReliablePurity ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        String mutationalLoadString = hasReliablePurity ? analysis().tumorMutationalLoadStatus().display() + " ("
                + SINGLE_DECIMAL_FORMAT.format(analysis().tumorMutationalBurden()) + " mut/genome)" : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Tumor mutational load").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(mutationalLoadString).addStyle(dataStyle)));

        String microSatelliteStabilityString = hasReliablePurity ? analysis().microsatelliteStatus().display() + " ("
                + DOUBLE_DECIMAL_FORMAT.format(analysis().microsatelliteIndelsPerMb()) + ")" : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Microsatellite (in)stability").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(microSatelliteStabilityString).addStyle(dataStyle)));

        String hrdString;
        Style hrdStyle;

        if (hasReliablePurity && (ChordStatus.HR_DEFICIENT == analysis().chordHrdStatus()
                || ChordStatus.HR_PROFICIENT == analysis().chordHrdStatus())) {
            hrdString = analysis().chordHrdStatus().display() + " (" + DOUBLE_DECIMAL_FORMAT.format(analysis().chordHrdValue()) + ")";
            hrdStyle = ReportResources.dataHighlightStyle();
        } else {
            hrdString = DataUtil.NA_STRING;
            hrdStyle = ReportResources.dataHighlightNaStyle();
        }

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("HR Status").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(hrdString).addStyle(hrdStyle)));

        Style virusStyle;
        String virusInterpretationString;
        String virusSummary = patientReport.sampleReport().reportViralInsertions()
                ? VirusBreakends.virusInterpretationSummary(analysis().virusBreakends())
                : DataUtil.NA_STRING;
        if (!virusSummary.equals(DataUtil.NA_STRING)) {
            virusStyle = ReportResources.dataHighlightStyle();
            virusInterpretationString = virusSummary;
        } else {
            virusStyle = ReportResources.dataHighlightNaStyle();
            virusInterpretationString = DataUtil.NA_STRING;
        }

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Integrated Virus").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(virusInterpretationString).addStyle(virusStyle)));
        div.add(table);

        reportDocument.add(div);
    }

    private static void renderTumorPurity(boolean hasReliablePurity, @NotNull String valueLabel, double value, double min, double max,
            @NotNull Table table) {
        String label = "Tumor purity";
        table.addCell(createMiddleAlignedCell().add(new Paragraph(label).addStyle(ReportResources.bodyTextStyle())));

        if (hasReliablePurity) {
            table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(valueLabel).addStyle(ReportResources.dataHighlightStyle())));
            table.addCell(createMiddleAlignedCell().add(createInlineBarChart(value, min, max)));
        } else {
            table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(Lims.PURITY_NOT_RELIABLE_STRING).addStyle(ReportResources.dataHighlightNaStyle())));
        }
    }

    private void renderGenomicAlterations(@NotNull Document report) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCellSummary()
                .add(new Paragraph("Genomic alterations in cancer genes").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        Set<String> driverVariantGenes = SomaticVariants.driverGenesWithVariant(analysis().reportableVariants());

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with driver mutation").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(driverVariantGenes)));

        int reportedVariants = SomaticVariants.countReportableVariants(analysis().reportableVariants());
        Style reportedVariantsStyle =
                (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of reported variants").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(String.valueOf(reportedVariants)).addStyle(
                reportedVariantsStyle)));

        Set<String> amplifiedGenes = GainsAndLosses.amplifiedGenes(analysis().gainsAndLosses());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Amplified gene(s)").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(amplifiedGenes)));

        Set<String> copyLossGenes = GainsAndLosses.lostGenes(analysis().gainsAndLosses());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Deleted gene(s)").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(copyLossGenes)));

        Set<String> disruptedGenes = HomozygousDisruptions.disruptedGenes(analysis().homozygousDisruptions());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Homozygously disrupted genes").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(disruptedGenes)));

        Set<String> fusionGenes = GeneFusions.uniqueGeneFusions(analysis().geneFusions());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Gene fusions").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(fusionGenes)));

        div.add(table);

        report.add(div);
    }

    private void renderPeach(@NotNull Document report) {
        Div div = createSectionStartDiv(contentWidth());
        String title = "Pharmacogenetics";

        Table table;
        if (patientReport.genomicAnalysis().hasReliablePurity() && patientReport.qsFormNumber().equals(QsFormNumber.FOR_080.display())) {

            table = TableUtil.createReportContentTableSummary(new float[] { 15, 20, 25, 50 },
                    new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Genotype"),
                            TableUtil.createHeaderCell("Function"),
                            TableUtil.createHeaderCell("Linked drugs").setTextAlignment(TextAlignment.CENTER) });

            for (PeachGenotype peachGenotype : Pharmacogenetics.sort(analysis().peachGenotypes())) {
                table.addCell(TableUtil.createContentCell(peachGenotype.gene()));
                table.addCell(TableUtil.createContentCell(peachGenotype.haplotype()));
                table.addCell(TableUtil.createContentCell(peachGenotype.function()));
                table.addCell(TableUtil.createContentCell(peachGenotype.linkedDrugs()).setTextAlignment(TextAlignment.CENTER));
            }
            table = TableUtil.createWrappingReportTableSummary(title, table);
        } else {
            table = TableUtil.createNAReportTable(title);
        }

        div.add(table);
        report.add(div);
    }

    @NotNull
    @VisibleForTesting
    static Set<String> sortGenes(@NotNull Set<String> driverVariantGenes) {
        List<String> genesList = Lists.newArrayList(driverVariantGenes);
        Collections.sort(genesList);
        return Sets.newHashSet(genesList);
    }

    @NotNull
    private static Div createSectionStartDiv(float width) {
        return new Div().setKeepTogether(true).setWidth(width).add(LineDivider.createLineDivider(width));
    }

    @NotNull
    private static Cell createMiddleAlignedCell() {
        return createMiddleAlignedCell(1);
    }

    @NotNull
    private static Cell createMiddleAlignedCell(int colSpan) {
        return TableUtil.createLayoutCell(1, colSpan).setVerticalAlignment(VerticalAlignment.MIDDLE);
    }

    @NotNull
    private static Cell createGeneListCell(@NotNull Set<String> genes) {
        String geneString = (genes.size() > 0) ? String.join(", ", genes) : DataUtil.NONE_STRING;

        Style style = (genes.size() > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        return createMiddleAlignedCell().add(createHighlightParagraph(geneString)).addStyle(style);
    }

    @NotNull
    private static Paragraph createHighlightParagraph(@NotNull String text) {
        return new Paragraph(text).setFixedLeading(14);
    }

    @NotNull
    private static Cell createTreatmentIndicationCell(int eventCount, int treatmentCount, @NotNull String treatmentsName) {
        String treatmentText;
        Style style;
        if (eventCount > 0) {
            treatmentText = String.format("%d | %d %s", eventCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = DataUtil.NONE_STRING;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell().add(createHighlightParagraph(treatmentText)).addStyle(style);
    }

    @NotNull
    private static InlineBarChart createInlineBarChart(double value, double min, double max) {
        InlineBarChart chart = new InlineBarChart(value, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }
}
