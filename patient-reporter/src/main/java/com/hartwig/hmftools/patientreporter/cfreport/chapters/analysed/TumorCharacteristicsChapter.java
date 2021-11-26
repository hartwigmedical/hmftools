package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import java.net.MalformedURLException;
import java.text.DecimalFormat;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.BarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.HrDeficiency;
import com.hartwig.hmftools.patientreporter.cfreport.data.MicroSatelliteStatus;
import com.hartwig.hmftools.patientreporter.cfreport.data.MutationalBurden;
import com.hartwig.hmftools.patientreporter.cfreport.data.MutationalLoad;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorCharacteristicsChapter implements ReportChapter {

    private static final float TABLE_SPACER_HEIGHT = 30;

    private static final DecimalFormat NO_DECIMAL_FORMAT = ReportResources.decimalFormat("#");
    private static final DecimalFormat SINGLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat DOUBLE_DECIMAL_FORMAT = ReportResources.decimalFormat("#.##");

    @NotNull
    private final AnalysedPatientReport patientReport;

    public TumorCharacteristicsChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    public String name() {
        return "Tumor characteristics";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        renderHrdCharacteristic(reportDocument);
        renderMicrosatelliteStabilityCharacteristic(reportDocument);
        renderMutationalLoadCharacteristic(reportDocument);
        renderMutationalBurdenCharacteristic(reportDocument);
        renderMolecularTissueOriginPlot(reportDocument);
    }

    private void renderHrdCharacteristic(@NotNull Document reportDocument) {
        GenomicAnalysis genomicAnalysis = patientReport.genomicAnalysis();
        double hrdValue = genomicAnalysis.chordHrdValue();
        ChordStatus hrdStatus = genomicAnalysis.chordHrdStatus();

        boolean hasReliablePurity = genomicAnalysis.hasReliablePurity();
        String hrDeficiencyLabel =
                hasReliablePurity ? hrdStatus.display() + " " + DOUBLE_DECIMAL_FORMAT.format(hrdValue) : DataUtil.NA_STRING;

        String hrdUnreliableFootnote = "* HRD score can not be determined reliably when a tumor is microsatellite unstable "
                + "(MSI) or has insufficient number of mutations and is therefore not reported for this sample.";
        boolean displayFootNote = false;
        boolean isHrdReliable = genomicAnalysis.chordHrdStatus() == ChordStatus.HR_PROFICIENT
                || genomicAnalysis.chordHrdStatus() == ChordStatus.HR_DEFICIENT;
        if (!isHrdReliable) {
            displayFootNote = true;
            hrDeficiencyLabel = DataUtil.NA_STRING + "*";
        }

        // We subtract 0.0001 from the minimum to allow visualization of a HR-score of exactly 0.
        BarChart hrChart = new BarChart(hrdValue, HrDeficiency.RANGE_MIN - 0.0001, HrDeficiency.RANGE_MAX, "Low", "High", false);
        hrChart.enabled(hasReliablePurity && isHrdReliable);
        hrChart.setTickMarks(HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, 0.1, SINGLE_DECIMAL_FORMAT);

        hrChart.setIndicator(ChordStatus.HRD_THRESHOLD, "HRD status (" + DOUBLE_DECIMAL_FORMAT.format(ChordStatus.HRD_THRESHOLD) + ")");

        reportDocument.add(createCharacteristicDiv("HR-Deficiency score",
                hrDeficiencyLabel,
                "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing "
                        + "the signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation. \n"
                        + "Tumors with a score greater or equal than 0.5 are considered HR deficient by complete BRCA inactivation.",
                hrChart,
                hrdUnreliableFootnote,
                displayFootNote));
    }

    private void renderMicrosatelliteStabilityCharacteristic(@NotNull Document reportDocument) {
        GenomicAnalysis genomicAnalysis = patientReport.genomicAnalysis();
        boolean hasReliablePurity = genomicAnalysis.hasReliablePurity();
        double microSatelliteStability = genomicAnalysis.microsatelliteIndelsPerMb();
        String microSatelliteStabilityString = hasReliablePurity ? genomicAnalysis.microsatelliteStatus().display() + " "
                + DOUBLE_DECIMAL_FORMAT.format(genomicAnalysis.microsatelliteIndelsPerMb()) : DataUtil.NA_STRING;

        BarChart satelliteChart =
                new BarChart(microSatelliteStability, MicroSatelliteStatus.RANGE_MIN, MicroSatelliteStatus.RANGE_MAX, "MSS", "MSI", false);
        satelliteChart.enabled(hasReliablePurity);
        satelliteChart.scale(InlineBarChart.LOG10_SCALE);
        satelliteChart.setTickMarks(new double[] { MicroSatelliteStatus.RANGE_MIN, 10, MicroSatelliteStatus.RANGE_MAX },
                DOUBLE_DECIMAL_FORMAT);
        satelliteChart.enableUndershoot(NO_DECIMAL_FORMAT.format(0));
        satelliteChart.enableOvershoot(">" + NO_DECIMAL_FORMAT.format(satelliteChart.max()));
        satelliteChart.setIndicator(MicroSatelliteStatus.THRESHOLD,
                "Microsatellite \ninstability (" + DOUBLE_DECIMAL_FORMAT.format(MicroSatelliteStatus.THRESHOLD) + ")");
        reportDocument.add(createCharacteristicDiv("Microsatellite status",
                microSatelliteStabilityString,
                "The microsatellite stability score represents the number of somatic inserts and deletes in "
                        + "(short) repeat sections across the whole genome of the tumor per Mb. This metric can be "
                        + "considered as a good marker for instability in microsatellite repeat regions. Tumors with a "
                        + "score greater than 4.0 are considered microsatellite unstable (MSI).",
                satelliteChart,
                Strings.EMPTY,
                false));
    }

    private void renderMutationalLoadCharacteristic(@NotNull Document reportDocument) {
        GenomicAnalysis genomicAnalysis = patientReport.genomicAnalysis();

        boolean hasReliablePurity = genomicAnalysis.hasReliablePurity();
        int mutationalLoad = genomicAnalysis.tumorMutationalLoad();
        String tmlStatus = genomicAnalysis.tumorMutationalLoadStatus().display();

        String mutationalLoadString = hasReliablePurity ? tmlStatus + " " + NO_DECIMAL_FORMAT.format(mutationalLoad) : DataUtil.NA_STRING;
        BarChart mutationalLoadChart =
                new BarChart(mutationalLoad, MutationalLoad.RANGE_MIN, MutationalLoad.RANGE_MAX, "Low", "High", false);
        mutationalLoadChart.enabled(hasReliablePurity);
        mutationalLoadChart.scale(InlineBarChart.LOG10_SCALE);
        mutationalLoadChart.setTickMarks(new double[] { MutationalLoad.RANGE_MIN, 10, 100, MutationalLoad.RANGE_MAX }, NO_DECIMAL_FORMAT);
        mutationalLoadChart.enableUndershoot(NO_DECIMAL_FORMAT.format(0));
        mutationalLoadChart.enableOvershoot(">" + NO_DECIMAL_FORMAT.format(mutationalLoadChart.max()));
        mutationalLoadChart.setIndicator(MutationalLoad.THRESHOLD,
                "Eligible for \nDRUP (" + NO_DECIMAL_FORMAT.format(MutationalLoad.THRESHOLD) + ")");

        reportDocument.add(createCharacteristicDiv("Tumor mutational load",
                mutationalLoadString,
                "The tumor mutational load represents the total number of somatic missense variants across "
                        + "the whole genome of the tumor. Patients with a mutational load over 140 could be eligible for "
                        + "immunotherapy within the DRUP study.",
                mutationalLoadChart,
                Strings.EMPTY,
                false));
    }

    private void renderMutationalBurdenCharacteristic(@NotNull Document reportDocument) {
        GenomicAnalysis genomicAnalysis = patientReport.genomicAnalysis();

        boolean hasReliablePurity = genomicAnalysis.hasReliablePurity();
        double mutationalBurden = genomicAnalysis.tumorMutationalBurden();
        String mutationalBurdenString =
                hasReliablePurity ? SINGLE_DECIMAL_FORMAT.format(mutationalBurden) + " variants per Mb" : DataUtil.NA_STRING;
        BarChart mutationalBurdenChart =
                new BarChart(mutationalBurden, MutationalBurden.RANGE_MIN, MutationalBurden.RANGE_MAX, "Low", "High", false);
        mutationalBurdenChart.enabled(hasReliablePurity);
        mutationalBurdenChart.scale(InlineBarChart.LOG10_SCALE);
        mutationalBurdenChart.setTickMarks(new double[] { MutationalBurden.RANGE_MIN, 10, MutationalBurden.RANGE_MAX },
                DOUBLE_DECIMAL_FORMAT);
        mutationalBurdenChart.enableUndershoot(NO_DECIMAL_FORMAT.format(0));
        mutationalBurdenChart.enableOvershoot(">" + SINGLE_DECIMAL_FORMAT.format(mutationalBurdenChart.max()));
        reportDocument.add(createCharacteristicDiv("Tumor mutational burden",
                mutationalBurdenString,
                "The tumor mutational burden score represents the number of all somatic variants across the "
                        + "whole genome of the tumor per Mb.",
                mutationalBurdenChart,
                Strings.EMPTY,
                false));
    }

    private void renderMolecularTissueOriginPlot(@NotNull Document reportDocument) {
        //TODO: fix how to add on new page?
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));
        reportDocument.add(createCharacteristicDiv(""));

        reportDocument.add(createCharacteristicDiv("Molecular tissue of origin prediction"));
        Table table = new Table(UnitValue.createPercentArray(new float[] { 10, 1, 10, 1, 10 }));
        table.setWidth(contentWidth());
        if (patientReport.molecularTissueOrigin() != null && patientReport.genomicAnalysis().hasReliablePurity()) {

            String molecularTissueOriginPlot = patientReport.molecularTissueOrigin().plotPath();
            if (patientReport.qsFormNumber().equals(QsFormNumber.FOR_209.display()) || patientReport.qsFormNumber()
                    .equals(QsFormNumber.FOR_080.display())) {
                if (patientReport.genomicAnalysis().impliedPurity() < 0.20) {
                    reportDocument.add(createCharacteristicDisclaimerDiv(
                            "Due to the low tumor purity, the molecular tissue of origin prediction should be interpreted with caution."));
                }

                try {
                    reportDocument.add(createCharacteristicDiv("")); // For better display plot
                    Image circosImage = new Image(ImageDataFactory.create(molecularTissueOriginPlot));
                    circosImage.setMaxHeight(250);
                    circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
                    circosImage.setMarginBottom(8);
                    reportDocument.add(circosImage);
                } catch (MalformedURLException e) {
                    throw new IOException("Failed to read molecular tissue origin plot image at " + molecularTissueOriginPlot);
                }
            }

            reportDocument.add(createCharacteristicDiv(""));
            reportDocument.add(createCharacteristicDiv(""));
            reportDocument.add(createCharacteristicDiv(""));
            reportDocument.add(createCharacteristicDiv(""));

            table.addCell(TableUtil.createLayoutCell()
                    .add(new Div().add(createContentParagraph("The title",
                            " shows the conclusion of the prediction of the molecular"
                                    + " tissue of origin. If none of the similarity predictions has a likelihood ≥80%, no reliable conclusion"
                                    + " can be drawn (‘results inconclusive’)."))));

            table.addCell(TableUtil.createLayoutCell());

            table.addCell(TableUtil.createLayoutCell()
                    .add(new Div().add(createContentParagraph("The left plot",
                            " shows the likelihoods (similarity) for all the origin "
                                    + "types analyzed by the molecular tissue of origin prediction tool. Only when the likelihood is ≥80% "
                                    + "(a peak in the green outer band of the plot), a reliable prediction (with >95% accuracy) can be drawn. "
                                    + "Lower likelihoods (<80%) suggest there is similarity with that tissue of origin, but this is less strong "
                                    + "and there is lower confidence."))));

            table.addCell(TableUtil.createLayoutCell());

            table.addCell(TableUtil.createLayoutCell()
                    .add(new Div().add(createContentParagraph("The right plot(s)",
                            " shows the breakdown of the strongest predicted "
                                    + "likelihood(s) into the contribution of the 1) SNV types (related to those used in Cosmic signatures), 2) "
                                    + "driver landscape and passenger characteristics (e.g. tumor-type specific drivers), and 3) somatic mutation "
                                    + "pattern (mutation distribution across the genome)."))));
        } else {
            reportDocument.add(new Paragraph("The molecular tissue of origin prediction is unreliable due to low tumor purity and "
                    + "therefore the results are not available.").addStyle(ReportResources.subTextStyle()));
        }

        reportDocument.add(table);
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String boldPart, @NotNull String regularPart) {
        return new Paragraph(boldPart).addStyle(ReportResources.subTextBoldStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(regularPart).addStyle(ReportResources.subTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private Div createCharacteristicDiv(@NotNull String title) {
        Div div = new Div();
        div.setKeepTogether(true);
        div.add(new Paragraph(title).addStyle(ReportResources.sectionTitleStyle()));

        return div;
    }

    @NotNull
    private Div createCharacteristicDisclaimerDiv(@NotNull String title) {
        Div div = new Div();
        div.setKeepTogether(true);

        div.add(new Paragraph(title).addStyle(ReportResources.smallBodyHeadingStyle()));
        return div;
    }

    @NotNull
    private Div createCharacteristicDiv(@NotNull String title, @NotNull String highlight, @NotNull String description,
            @NotNull BarChart chart, @NotNull String footnote, boolean displayFootnote) {
        Div div = new Div();
        div.setKeepTogether(true);

        div.add(new Paragraph(title).addStyle(ReportResources.sectionTitleStyle()));

        Table table = new Table(UnitValue.createPercentArray(new float[] { 10, 1, 19 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(DataLabel.createDataLabel(highlight)));

        table.addCell(TableUtil.createLayoutCell(2, 1));

        table.addCell(TableUtil.createLayoutCell(2, 1).add(chart));
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph(description).addStyle(ReportResources.bodyTextStyle())
                        .setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableUtil.createLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT));
        div.add(table);

        if (displayFootnote) {
            div.add(new Paragraph(footnote).addStyle(ReportResources.subTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING));
        }

        return div;
    }
}
