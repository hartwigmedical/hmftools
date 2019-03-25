package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
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
    protected void renderChapterContent(Document report) {

        final DecimalFormat wholeNumberFormat = new DecimalFormat("#");
        final DecimalFormat oneDecimalFormat = new DecimalFormat("#.#");

        final float hrDeficiency = 0f;

        BarChart hrChart = new BarChart(hrDeficiency, 0.0f, 1.0f, "Low", "High", BarChart.createTickMarkLabels(0f, 1f, 0.1f, oneDecimalFormat));
        hrChart.setIndicator(0.5f, "HR-Deficient");
        report.add(createCharacteristicDiv("HR-Deficiency score", String.format("%.0f", hrDeficiency), "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing the signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation.", hrChart));

        final float microSatelliteStability = 1.5f;
        final String microSatelliteStabilityString = "Stable";
        BarChart satelliteChart = new BarChart(microSatelliteStability, 0f, 10f, "MSS", "MSI", BarChart.createTickMarkLabels(0f, 10f, 1f, wholeNumberFormat));
        satelliteChart.enableDefaultRangeOvershoot(">" + wholeNumberFormat.format(satelliteChart.getMax()));
        satelliteChart.setIndicator(4f, "Microsattelite unstable");
        report.add(createCharacteristicDiv("Microsatellite status", String.format("%s %.4f", microSatelliteStabilityString, microSatelliteStability), "The microsatellite stability score represents the number of somatic inserts and deletes in (short) repeat sections across the whole genome of the tumor per Mb. This metric can be considered as a good marker for instability in microsatellite repeat regions. Tumors with a score greater than 4.0 are considered microsatellite unstable (MSI).", satelliteChart));

        final float mutationalLoad = 250; //182;
        final String mutationalLoadString = "High";
        BarChart mutationalLoadChart = new BarChart(mutationalLoad, 0, 250, "Low", "High", BarChart.createTickMarkLabels(0f, 250, 25f, wholeNumberFormat));
        mutationalLoadChart.enableDefaultRangeOvershoot(">" + wholeNumberFormat.format(mutationalLoadChart.getMax()));
        mutationalLoadChart.setIndicator(140, "Eligible for DRUP");
        report.add(createCharacteristicDiv("Tumor mutational load", String.format("%s %.0f", mutationalLoadString, mutationalLoad), "The tumor mutational load represents the total number of somatic missense variants across the whole genome of the tumor. Patients with a mutational load over 140 could be eligible for immunotherapy within the DRUP study.", mutationalLoadChart));

        final float mutationalBurden = 30f; //13.6f;
        BarChart mutationalBurdenChart = new BarChart(mutationalBurden, 0, 25f, "Low", "High", BarChart.createTickMarkLabels(0f, 25, 2.5f, oneDecimalFormat));
        mutationalBurdenChart.enableDefaultRangeOvershoot(">" + oneDecimalFormat.format(mutationalBurdenChart.getMax()));
        report.add(createCharacteristicDiv("Tumor mutational burden", String.format("%.1f variants per Mb", mutationalBurden), "The tumor mutational burden score represents the number of all somatic variants across the whole genome of the tumor per Mb.", mutationalBurdenChart));

    }

    @NotNull
    private final Div createCharacteristicDiv(@NotNull String title, @NotNull String highlight, @NotNull String description, @NotNull BarChart chart) {

        // Initialize div
        Div div = new Div();
        div.setKeepTogether(true);

        // Add title
        div.add(new Paragraph(title)
                .addStyle(ReportResources.sectionTitleStyle()));

        // Add content table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 2}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell().add(DataLabel.createDataLabel(highlight)));
        table.addCell(TableHelper.getLayoutCell(2, 1).add(chart));
        table.addCell(TableHelper.getLayoutCell().add(new Paragraph(description).addStyle(ReportResources.bodyTextStyle())));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer
        div.add(table);

        return div;

    }



}
