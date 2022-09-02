package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.hartwig.hmftools.patientreporter.cfreport.data.GainsAndLosses;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.HLAAllele;
import com.hartwig.hmftools.patientreporter.cfreport.data.HomozygousDisruptions;
import com.hartwig.hmftools.patientreporter.cfreport.data.Pharmacogenetics;
import com.hartwig.hmftools.patientreporter.cfreport.data.SomaticVariants;
import com.hartwig.hmftools.patientreporter.cfreport.data.TumorPurity;
import com.hartwig.hmftools.patientreporter.cfreport.data.ViralPresence;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SummaryChapter implements ReportChapter {
    private static final Logger LOGGER = LogManager.getLogger(SummaryChapter.class);

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
    public String pdfTitle() {
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

    @NotNull
    @Override
    public String name() {
        return "Summary";
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
        reportDocument.add(TumorLocationAndTypeTable.createBiopsyLocationAndTumorLocation(patientReport.sampleReport()
                .primaryTumorLocationString(), patientReport.sampleReport().biopsyLocationString(), contentWidth()));
        reportDocument.add(new Paragraph());
        reportDocument.add(TumorLocationAndTypeTable.createTumorType(patientReport.sampleReport().primaryTumorTypeString(),
                contentWidth()));
        reportDocument.add(new Paragraph("\nThe information regarding 'primary tumor location', 'primary tumor type' and 'biopsy location'"
                + "  \nis based on information received from the originating hospital.").addStyle(ReportResources.subTextStyle()));

        renderClinicalConclusionText(reportDocument);
        renderSpecialRemarkText(reportDocument);
        renderTreatmentIndications(reportDocument);
        renderTumorCharacteristics(reportDocument);
        renderGenomicAlterations(reportDocument);
        renderPeach(reportDocument);
        renderHla(reportDocument);
    }

    private void renderClinicalConclusionText(@NotNull Document reportDocument) {
        String text = patientReport.clinicalSummary();
        String clinicalConclusion = Strings.EMPTY;
        if (text == null) {
            String sentence = "An overview of all detected oncogenic DNA aberrations can be found in the report";

            if (!analysis().hasReliablePurity()) {
                clinicalConclusion = "Of note, WGS analysis indicated a very low abundance of genomic aberrations, which can be caused "
                        + "by a low tumor percentage in the received tumor material or due to genomic very stable/normal tumor type. "
                        + "As a consequence no reliable tumor purity assessment is possible and no information regarding "
                        + "mutation copy number and tVAF can be provided.\n" + sentence;
            } else if (analysis().impliedPurity() < ReportResources.PURITY_CUTOFF) {
                double impliedPurityPercentage =
                        MathUtil.mapPercentage(analysis().impliedPurity(), TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
                clinicalConclusion = "Due to the lower sensitivity (" + DataUtil.formatPercentage(impliedPurityPercentage) + ") "
                        + "of this test potential (subclonal) DNA aberrations might not have been detected using this test. " + ""
                        + "This result should therefore be considered with caution.\n" + sentence;
            }
        } else {
            clinicalConclusion = text;
        }

        if (!clinicalConclusion.isEmpty()) {
            Div div = createSectionStartDiv(contentWidth());
            div.add(new Paragraph("Summary of clinical relevance").addStyle(ReportResources.sectionTitleStyle()));

            div.add(new Paragraph(text).setWidth(contentWidth()).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(11));
            div.add(new Paragraph("\nThis summary is generated based on DNA analysis only. Clinical patient characteristics or "
                    + "medical records have not been considered.").addStyle(ReportResources.subTextStyle()));

            reportDocument.add(div);
        }
    }

    private void renderSpecialRemarkText(@NotNull Document reportDocument) {
        String text = patientReport.specialRemark();

        if (!text.isEmpty()) {
            Div div = createSectionStartDiv(contentWidth());
            div.add(new Paragraph("Special Remark").addStyle(ReportResources.sectionTitleStyle()));

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
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with therapy indication").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(therapyEventCount,
                EvidenceItems.onLabelTreatmentString(analysis().tumorSpecificEvidence()),
                "treatment(s)"));

        int trialEventCount = ClinicalTrials.uniqueEventCount(analysis().clinicalTrials());
        int trialCount = ClinicalTrials.uniqueTrialCount(analysis().clinicalTrials());
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with clinical trial eligibility").addStyle(
                ReportResources.bodyTextStyle())));
        table.addCell(createStudyIndicationCell(trialEventCount, trialCount, "trial(s)"));
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

        String cuppaPrediction = patientReport.cuppaPrediction() != null && patientReport.genomicAnalysis().hasReliablePurity()
                ? patientReport.cuppaPrediction().cancerType() + " (" + patientReport.cuppaPrediction().likelihood() + ")"
                : DataUtil.NA_STRING;
        Style dataStyleMolecularTissuePrediction =
                hasReliablePurity ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Molecular tissue of origin prediction").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(cuppaPrediction).addStyle(dataStyleMolecularTissuePrediction)));

        Style dataStyle = hasReliablePurity ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        String mutationalLoadString = hasReliablePurity ? analysis().tumorMutationalLoadStatus().display() + " ("
                + SINGLE_DECIMAL_FORMAT.format(analysis().tumorMutationalLoad()) + " mut/genome)" : DataUtil.NA_STRING;
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

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Integrated Virus").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createVirusInterpretationString(ViralPresence.virusInterpretationSummary(analysis().reportableViruses()),
                patientReport.sampleReport().reportViralPresence()));

        div.add(table);

        reportDocument.add(div);
    }

    @NotNull
    private static Cell createVirusInterpretationString(@NotNull Set<String> virus, boolean reportViralPresence) {
        String virusSummary;
        Style style;
        if (reportViralPresence && virus.size() == 0) {
            virusSummary = DataUtil.NONE_STRING;
            style = ReportResources.dataHighlightNaStyle();
        } else if (reportViralPresence && virus.size() > 0) {
            virusSummary = String.join(", ", virus);
            style = ReportResources.dataHighlightStyle();
        } else {
            virusSummary = DataUtil.NA_STRING;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell(2).add(createHighlightParagraph(virusSummary)).addStyle(style);
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

        MicrosatelliteStatus microSatelliteStabilityString =
                analysis().hasReliablePurity() ? analysis().microsatelliteStatus() : MicrosatelliteStatus.UNKNOWN;
        if (microSatelliteStabilityString == MicrosatelliteStatus.MSI) {
            Set<String> genesDisplay = SomaticVariants.determineMSIgenes(analysis().reportableVariants(),
                    analysis().gainsAndLosses(),
                    analysis().homozygousDisruptions());
            table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                    .add(new Paragraph("Potential MMR genes").addStyle(ReportResources.bodyTextStyle())));
            table.addCell(createGeneListCell(sortGenes(genesDisplay)));
        }

        ChordStatus chordStatus = analysis().hasReliablePurity() ? analysis().chordHrdStatus() : ChordStatus.UNKNOWN;
        if (chordStatus == ChordStatus.HR_DEFICIENT) {
            Set<String> genesDisplay = SomaticVariants.determineHRDgenes(analysis().reportableVariants(),
                    analysis().gainsAndLosses(),
                    analysis().homozygousDisruptions());
            table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                    .add(new Paragraph("Potential HRD genes").addStyle(ReportResources.bodyTextStyle())));
            table.addCell(createGeneListCell(sortGenes(genesDisplay)));
        }

        div.add(table);

        report.add(div);
    }

    private void renderPeach(@NotNull Document report) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCellSummary()
                .add(new Paragraph("Pharmacogenetics").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        Set<String> pgxFunctions;
        Set<String> pgxGenes;
        Style pgxStyle;
        String reportedPhenotypes;

        if (patientReport.sampleReport().reportPharmogenetics() && patientReport.peachGenotypes().size() > 0) {
            pgxFunctions = Pharmacogenetics.phenotypesFunctions(patientReport.peachGenotypes());
            pgxGenes = Pharmacogenetics.phenotypesGenes(patientReport.peachGenotypes());
            pgxStyle = ReportResources.dataHighlightStyle();
            reportedPhenotypes = Integer.toString(Pharmacogenetics.countPhenotypes(patientReport.peachGenotypes()));
        } else {
            pgxFunctions = Sets.newHashSet(DataUtil.NA_STRING);
            pgxGenes = Sets.newHashSet(DataUtil.NA_STRING);
            pgxStyle = ReportResources.dataHighlightNaStyle();
            reportedPhenotypes = DataUtil.NA_STRING;
        }

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with haplotypes").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(pgxGenes)).addStyle(pgxStyle));

        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of reported haplotypes").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(reportedPhenotypes).addStyle(pgxStyle)));

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Functions of the haplotypes").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(pgxFunctions)).addStyle(pgxStyle));

        div.add(table);

        report.add(div);
    }

    private void renderHla(@NotNull Document report) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCellSummary().add(new Paragraph("HLA Alleles").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        Set<String> HLAtypes = Sets.newHashSet();
        Set<String> HLBtypes = Sets.newHashSet();
        Set<String> HLCtypes = Sets.newHashSet();
        Style hlaStyle = ReportResources.dataHighlightStyle();

        for (LilacAllele lilacAllele : patientReport.genomicAnalysis().lilac().alleles()) {
            String allele = lilacAllele.allele();
            if (allele.startsWith("A*")) {
                HLAtypes.add(allele.substring(1, allele.length()-1));
            }
            else if (allele.startsWith("B*")) {
                HLBtypes.add(allele.substring(1, allele.length()-1));
            }
            else if (allele.startsWith("C*")) {
                HLCtypes.add(allele.substring(1, allele.length()-1));
            }
        }

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("HLA Alleles").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(HLAtypes).addStyle(hlaStyle));

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("HLB Alleles").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(HLBtypes).addStyle(hlaStyle));

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("HLC Alleles").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(HLCtypes).addStyle(hlaStyle));

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
    private static Cell createTreatmentIndicationCell(int eventCount, @NotNull String treatmentCount, @NotNull String treatmentsName) {
        String treatmentText;
        Style style;
        if (eventCount > 0) {
            treatmentText = String.format("%d | %s %s", eventCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = DataUtil.NONE_STRING;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell().add(createHighlightParagraph(treatmentText)).addStyle(style);
    }

    @NotNull
    private static Cell createStudyIndicationCell(int eventCount, int treatmentCount, @NotNull String treatmentsName) {
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
