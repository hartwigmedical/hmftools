package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.*;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

public class TumorCharacteristicsChapter implements ReportChapter {

    private final static float TABLE_SPACER_HEIGHT = 30;

    private final AnalysedPatientReport patientReport;

    public TumorCharacteristicsChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String getName() {
        return "Tumor characteristics";
    }

    @Override
    public final void render(@NotNull final Document reportDocument) {

        final DecimalFormat noDecimalFormat = new DecimalFormat("#");
        final DecimalFormat singleDecimalFormat = new DecimalFormat("#.#");
        final DecimalFormat doubleDecimalFormat = new DecimalFormat("#.##");

        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        // HR Deficiency
        final double hrDeficiency = patientReport.chordAnalysis().hrdValue();
        final String hrDeficiencyLabel = HrDeficiency.interpretToString(hrDeficiency, hasReliablePurityFit);
        BarChart hrChart = new BarChart(hrDeficiency, HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, "Low", "High");
        hrChart.setEnabled(hasReliablePurityFit);
        hrChart.setTickMarks(HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, 0.1, singleDecimalFormat);
        // @TODO HR threshold to be determined: hrChart.setIndicator(0.5f, "HR-Deficient");
        reportDocument.add(createCharacteristicDiv("HR-Deficiency score", hrDeficiencyLabel,
                "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing " +
                        "the signature of this sample with signatures found across samples with known BRCA1/BRCA2 " +
                        "inactivation.", hrChart));


        // Microsatellite stability
        final double microSatelliteStability = patientReport.microsatelliteIndelsPerMb();
        final String microSatelliteStabilityString = hasReliablePurityFit
                ? MicroSatelliteStatus.interpretToString(microSatelliteStability) + " " + new DecimalFormat("#.####").format(microSatelliteStability)
                : DataUtil.NAString;
        BarChart satelliteChart = new BarChart(
                microSatelliteStability,
                MicroSatelliteStatus.RANGE_MIN,
                MicroSatelliteStatus.RANGE_MAX,
                "MSS", "MSI");
        satelliteChart.setEnabled(hasReliablePurityFit);
        satelliteChart.setScale(InlineBarChart.LOG10_SCALE);
        satelliteChart.setTickMarks(new double[] {MicroSatelliteStatus.RANGE_MIN, 10, MicroSatelliteStatus.RANGE_MAX}, doubleDecimalFormat);
        satelliteChart.enableUndershoot("<" + noDecimalFormat.format(satelliteChart.getMin()));
        satelliteChart.enableOvershoot(">" + noDecimalFormat.format(satelliteChart.getMax()));
        satelliteChart.setIndicator(MicroSatelliteStatus.THRESHOLD, "Microsatellite \ninstability (" + doubleDecimalFormat.format(MicroSatelliteStatus.THRESHOLD) + ")");
        reportDocument.add(createCharacteristicDiv("Microsatellite status", microSatelliteStabilityString,
                "The microsatellite stability score represents the number of somatic inserts and deletes in " +
                        "(short) repeat sections across the whole genome of the tumor per Mb. This metric can be " +
                        "considered as a good marker for instability in microsatellite repeat regions. Tumors with a " +
                        "score greater than 4.0 are considered microsatellite unstable (MSI).", satelliteChart));


        // Mutational load
        final int mutationalLoad = patientReport.tumorMutationalLoad();
        final String mutationalLoadString = hasReliablePurityFit
                ? MutationalLoad.interpretToString(mutationalLoad, hasReliablePurityFit) + " " + noDecimalFormat.format(mutationalLoad)
                : DataUtil.NAString;
        BarChart mutationalLoadChart = new BarChart(
                mutationalLoad,
                MutationalLoad.RANGE_MIN,
                MutationalLoad.RANGE_MAX,
                "Low", "High");
        mutationalLoadChart.setEnabled(hasReliablePurityFit);
        mutationalLoadChart.setScale(InlineBarChart.LOG10_SCALE);
        mutationalLoadChart.setTickMarks(
                new double[] {MutationalLoad.RANGE_MIN, 10, 100, MutationalLoad.RANGE_MAX}, noDecimalFormat);
        mutationalLoadChart.enableUndershoot("<" + noDecimalFormat.format(mutationalLoadChart.getMin()));
        mutationalLoadChart.enableOvershoot(">" + noDecimalFormat.format(mutationalLoadChart.getMax()));
        mutationalLoadChart.setIndicator(MutationalLoad.THRESHOLD, "Eligible for \nDRUP (" + noDecimalFormat.format(MutationalLoad.THRESHOLD) + ")");

        reportDocument.add(createCharacteristicDiv("Tumor mutational load", mutationalLoadString,
                "The tumor mutational load represents the total number of somatic missense variants across " +
                        "the whole genome of the tumor. Patients with a mutational load over 140 could be eligible for " +
                        "immunotherapy within the DRUP study.", mutationalLoadChart));


        // Mutational burden
        final double mutationalBurden = patientReport.tumorMutationalBurden();
        final String mutationalBurdenString = hasReliablePurityFit
                ? singleDecimalFormat.format(mutationalBurden) + " variants per Mb"
                : DataUtil.NAString;
        BarChart mutationalBurdenChart = new BarChart(
                mutationalBurden,
                MutationalBurden.RANGE_MIN,
                MutationalBurden.RANGE_MAX,
                "Low", "High");
        mutationalBurdenChart.setEnabled(hasReliablePurityFit);
        mutationalBurdenChart.setScale(InlineBarChart.LOG10_SCALE);
        mutationalBurdenChart.setTickMarks(new double[] {MutationalBurden.RANGE_MIN, 10, MutationalBurden.RANGE_MAX}, doubleDecimalFormat);
        mutationalBurdenChart.enableUndershoot("<" + singleDecimalFormat.format(mutationalBurdenChart.getMin()));
        mutationalBurdenChart.enableOvershoot(">" + singleDecimalFormat.format(mutationalBurdenChart.getMax()));
        reportDocument.add(createCharacteristicDiv("Tumor mutational burden", mutationalBurdenString,
                "The tumor mutational burden score represents the number of all somatic variants across the " +
                        "whole genome of the tumor per Mb.", mutationalBurdenChart));


    }

    @NotNull
    private Div createCharacteristicDiv(@NotNull String title, @NotNull String highlight, @NotNull String description, @NotNull BarChart chart) {

        // Initialize div
        Div div = new Div();
        div.setKeepTogether(true);

        // Add title
        div.add(new Paragraph(title)
                .addStyle(ReportResources.sectionTitleStyle()));

        // Add content table
        Table table = new Table(UnitValue.createPercentArray(new float[] {10, 1, 19}));
        table.setWidth(getContentWidth());
        table.addCell(TableUtil.getLayoutCell().add(DataLabel.createDataLabel(highlight)));

        table.addCell(TableUtil.getLayoutCell(2,1)); // Spacer

        table.addCell(TableUtil.getLayoutCell(2, 1).add(chart));
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph(description)
                        .addStyle(ReportResources.bodyTextStyle())
                        .setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableUtil.getLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT)); // Spacer
        div.add(table);

        return div;

    }

}
