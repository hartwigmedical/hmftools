package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.DataUtility;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;

public class TumorCharacteristicsChapter extends ReportChapter {

    private final static float TABLE_SPACER_HEIGHT = 30;

    @Override
    public String getName() {
        return "Tumor characteristics";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {

        final DecimalFormat noDecimalFormat = new DecimalFormat("#");
        final DecimalFormat singleDecimalFormat = new DecimalFormat("#.#");
        final DecimalFormat doubleDecimalFormat = new DecimalFormat("#.##");


        // HR Deficiency
        final double hrDeficiency = patientReport.chordAnalysis().hrdValue();
        final String hrDeficiencyLabel = patientReport.hasReliablePurityFit()
                ? DataUtility.HrDeficiency.interpretToString(hrDeficiency)
                : DataUtility.NAString;
        BarChart hrChart = new BarChart(hrDeficiency, 0.0f, 1.0f, "Low", "High");
        hrChart.setTickMarks(0, 1, 0.1, singleDecimalFormat);
        hrChart.setIndicator(0.5f, "HR-Deficient");
        reportDocument.add(createCharacteristicDiv("HR-Deficiency score", hrDeficiencyLabel,
                "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing " +
                        "the signature of this sample with signatures found across samples with known BRCA1/BRCA2 " +
                        "inactivation.", hrChart));


        // Microsatellite stability
        final double microSatelliteStability = patientReport.microsatelliteIndelsPerMb();
        final String microSatelliteStabilityString = patientReport.hasReliablePurityFit()
                ? DataUtility.MicroSatellite.interpretToString(microSatelliteStability) + " " + new DecimalFormat("#.####").format(microSatelliteStability)
                : DataUtility.NAString;
        BarChart satelliteChart = new BarChart(
                microSatelliteStability,
                DataUtility.MicroSatellite.RANGE_MIN,
                DataUtility.MicroSatellite.RANGE_MAX,
                "MSS", "MSI");
        satelliteChart.setScale(InlineBarChart.LOG10_SCALE);
        satelliteChart.setTickMarks(new double[] {DataUtility.MicroSatellite.RANGE_MIN, 1E-1, 1, 10, DataUtility.MicroSatellite.RANGE_MAX}, doubleDecimalFormat);
        satelliteChart.enableDefaultRangeOvershoot(">" + noDecimalFormat.format(satelliteChart.getMax()));
        satelliteChart.setIndicator(DataUtility.MicroSatellite.THRESHOLD, "Microsatellite unstable");
        reportDocument.add(createCharacteristicDiv("Microsatellite status", microSatelliteStabilityString,
                "The microsatellite stability score represents the number of somatic inserts and deletes in " +
                        "(short) repeat sections across the whole genome of the tumor per Mb. This metric can be " +
                        "considered as a good marker for instability in microsatellite repeat regions. Tumors with a " +
                        "score greater than 4.0 are considered microsatellite unstable (MSI).", satelliteChart));


        // Mutational load
        final int mutationalLoad = patientReport.tumorMutationalLoad();
        final String mutationalLoadString = patientReport.hasReliablePurityFit()
                ? DataUtility.MutationalLoad.interpretToString(mutationalLoad) + " " + noDecimalFormat.format(mutationalLoad)
                : DataUtility.NAString;
        BarChart mutationalLoadChart = new BarChart(
                mutationalLoad,
                DataUtility.MutationalLoad.RANGE_MIN,
                DataUtility.MutationalLoad.RANGE_MAX,
                "Low", "High");
        mutationalLoadChart.setScale(InlineBarChart.LOG10_SCALE);
        mutationalLoadChart.setTickMarks(
                new double[] {DataUtility.MutationalLoad.RANGE_MIN, 10, 100, DataUtility.MutationalLoad.RANGE_MAX}, noDecimalFormat);
        mutationalLoadChart.enableDefaultRangeOvershoot(">" + noDecimalFormat.format(mutationalLoadChart.getMax()));
        mutationalLoadChart.setIndicator(DataUtility.MutationalLoad.THRESHOLD, "Eligible for DRUP");
        reportDocument.add(createCharacteristicDiv("Tumor mutational load", mutationalLoadString,
                "The tumor mutational load represents the total number of somatic missense variants across " +
                        "the whole genome of the tumor. Patients with a mutational load over 140 could be eligible for " +
                        "immunotherapy within the DRUP study.", mutationalLoadChart));



        // Mutational burden
        final double mutationalBurden = patientReport.tumorMutationalBurden();
        final String mutationalBurdenString = patientReport.hasReliablePurityFit()
                ? singleDecimalFormat.format(mutationalBurden) + " variants per Mb"
                : DataUtility.NAString;
        BarChart mutationalBurdenChart = new BarChart(
                mutationalBurden,
                DataUtility.MutationalBurden.RANGE_MIN,
                DataUtility.MutationalBurden.RANGE_MAX,
                "Low", "High");
        mutationalBurdenChart.setScale(InlineBarChart.LOG10_SCALE);
        mutationalBurdenChart.setTickMarks(new double[] {DataUtility.MutationalBurden.RANGE_MIN, 1E-1, 1, 10, DataUtility.MutationalBurden.RANGE_MAX}, doubleDecimalFormat);
        mutationalBurdenChart.enableDefaultRangeOvershoot(">" + singleDecimalFormat.format(mutationalBurdenChart.getMax()));
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
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 0.1f, 1.9f}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell().add(DataLabel.createDataLabel(highlight)));

        table.addCell(TableHelper.getLayoutCell(2,1)); // Spacer

        table.addCell(TableHelper.getLayoutCell(2, 1).add(chart));
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph(description)
                        .addStyle(ReportResources.bodyTextStyle())
                        .setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer
        div.add(table);

        return div;

    }



}
