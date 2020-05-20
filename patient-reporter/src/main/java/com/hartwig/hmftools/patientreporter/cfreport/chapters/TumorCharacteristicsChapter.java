package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.text.DecimalFormat;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.HrDeficiency;
import com.hartwig.hmftools.patientreporter.cfreport.data.MicroSatelliteStatus;
import com.hartwig.hmftools.patientreporter.cfreport.data.MutationalBurden;
import com.hartwig.hmftools.patientreporter.cfreport.data.MutationalLoad;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

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
    }

    private void renderHrdCharacteristic(@NotNull Document reportDocument) {
        boolean hasReliablePurity = patientReport.hasReliablePurity();
        boolean isMicrosatelliteStable = patientReport.microsatelliteStatus() == MicrosatelliteStatus.MSS;

        double hrdValue = patientReport.chordHrdValue();
        ChordStatus hrdStatus = patientReport.chordHrdStatus();
        String hrDeficiencyLabel =
                hasReliablePurity ? hrdStatus.display() + " " + DOUBLE_DECIMAL_FORMAT.format(hrdValue) : DataUtil.NA_STRING;

        String noHrdWhenMicrosatelliteUnstable = "* HRD score can not be determined reliably when a tumor is microsatellite unstable (MSI) "
                + "and is therefore not reported for this sample.";
        boolean displayFootNote = false;
        if (!isMicrosatelliteStable) {
            displayFootNote = true;
            hrDeficiencyLabel = DataUtil.NA_STRING + "*";
        }

        BarChart hrChart = new BarChart(hrdValue, HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, "Low", "High", false);
        if (hrdValue == 0.0) {
            hrChart = new BarChart(hrdValue, HrDeficiency.RANGE_MIN_BARPLOT, HrDeficiency.RANGE_MAX, "Low", "High", false);
            hrChart.enabled(hasReliablePurity && isMicrosatelliteStable);
            hrChart.setTickMarks(HrDeficiency.RANGE_MIN_BARPLOT, HrDeficiency.RANGE_MAX, 0.1, SINGLE_DECIMAL_FORMAT);
            hrChart.enableUndershoot(NO_DECIMAL_FORMAT.format(0));
        } else {
            hrChart.enabled(hasReliablePurity && isMicrosatelliteStable);
            hrChart.setTickMarks(HrDeficiency.RANGE_MIN_BARPLOT_NO_UNDERSHOOT, HrDeficiency.RANGE_MAX, 0.1, SINGLE_DECIMAL_FORMAT);
        }

        hrChart.setIndicator(ChordStatus.HRD_THRESHOLD, "HRD status (" + DOUBLE_DECIMAL_FORMAT.format(ChordStatus.HRD_THRESHOLD) + ")");

        reportDocument.add(createCharacteristicDiv("HR-Deficiency score",
                hrDeficiencyLabel,
                "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing "
                        + "the signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation. \n"
                        + "Tumors with a score greater or equals than 0.5 are considered HR deficient by complete BRCA inactivation.",
                hrChart,
                noHrdWhenMicrosatelliteUnstable,
                displayFootNote));
    }

    private void renderMicrosatelliteStabilityCharacteristic(@NotNull Document reportDocument) {
        boolean hasReliablePurity = patientReport.hasReliablePurity();
        double microSatelliteStability = patientReport.microsatelliteIndelsPerMb();
        String microSatelliteStabilityString = hasReliablePurity ? patientReport.microsatelliteStatus().display() + " "
                + DOUBLE_DECIMAL_FORMAT.format(patientReport.microsatelliteIndelsPerMb()) : DataUtil.NA_STRING;

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
        boolean hasReliablePurity = patientReport.hasReliablePurity();
        int mutationalLoad = patientReport.tumorMutationalLoad();
        TumorMutationalStatus tmlStatus = patientReport.tumorMutationalLoadStatus();

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
        boolean hasReliablePurity = patientReport.hasReliablePurity();
        double mutationalBurden = patientReport.tumorMutationalBurden();
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
