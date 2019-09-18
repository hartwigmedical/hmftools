package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.text.DecimalFormat;

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

import org.jetbrains.annotations.NotNull;

public class TumorCharacteristicsChapter implements ReportChapter {

    private final static float TABLE_SPACER_HEIGHT = 30;

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
    public final void render(@NotNull final Document reportDocument) {
        DecimalFormat noDecimalFormat = ReportResources.decimalFormat("#");
        DecimalFormat singleDecimalFormat = ReportResources.decimalFormat("#.#");
        DecimalFormat doubleDecimalFormat = ReportResources.decimalFormat("#.##");

        boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        double hrDeficiency = patientReport.chordAnalysis().hrdValue();
        String hrDeficiencyLabel = hasReliablePurityFit ? HrDeficiency.interpretToString(hrDeficiency) : DataUtil.NA_STRING;
        BarChart hrChart = new BarChart(hrDeficiency, HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, "Low", "High", hasReliablePurityFit);
        hrChart.enabled(hasReliablePurityFit);
        hrChart.setTickMarks(HrDeficiency.RANGE_MIN, HrDeficiency.RANGE_MAX, 0.1, singleDecimalFormat);
        reportDocument.add(createCharacteristicDiv("HR-Deficiency score",
                hrDeficiencyLabel,
                "The HR-deficiency score is determined by CHORD, a WGS signature-based classifier comparing "
                        + "the signature of this sample with signatures found across samples with known BRCA1/BRCA2 inactivation.",
                hrChart));

        double microSatelliteStability = patientReport.microsatelliteIndelsPerMb();

        String microSatelliteStabilityString =
                hasReliablePurityFit ? MicroSatelliteStatus.interpretToString(microSatelliteStability) + " " + doubleDecimalFormat.format(
                        microSatelliteStability) : DataUtil.NA_STRING;
        BarChart satelliteChart =
                new BarChart(microSatelliteStability, MicroSatelliteStatus.RANGE_MIN, MicroSatelliteStatus.RANGE_MAX, "MSS", "MSI", false);
        satelliteChart.enabled(hasReliablePurityFit);
        satelliteChart.scale(InlineBarChart.LOG10_SCALE);
        satelliteChart.setTickMarks(new double[] { MicroSatelliteStatus.RANGE_MIN, 10, MicroSatelliteStatus.RANGE_MAX },
                doubleDecimalFormat);
        satelliteChart.enableUndershoot(noDecimalFormat.format(0));
        satelliteChart.enableOvershoot(">" + noDecimalFormat.format(satelliteChart.max()));
        satelliteChart.setIndicator(MicroSatelliteStatus.THRESHOLD,
                "Microsatellite \ninstability (" + doubleDecimalFormat.format(MicroSatelliteStatus.THRESHOLD) + ")");
        reportDocument.add(createCharacteristicDiv("Microsatellite status",
                microSatelliteStabilityString,
                "The microsatellite stability score represents the number of somatic inserts and deletes in "
                        + "(short) repeat sections across the whole genome of the tumor per Mb. This metric can be "
                        + "considered as a good marker for instability in microsatellite repeat regions. Tumors with a "
                        + "score greater than 4.0 are considered microsatellite unstable (MSI).",
                satelliteChart));

        int mutationalLoad = patientReport.tumorMutationalLoad();
        String mutationalLoadString = hasReliablePurityFit
                ? MutationalLoad.interpretToString(mutationalLoad) + " " + noDecimalFormat.format(mutationalLoad)
                : DataUtil.NA_STRING;
        BarChart mutationalLoadChart =
                new BarChart(mutationalLoad, MutationalLoad.RANGE_MIN, MutationalLoad.RANGE_MAX, "Low", "High", false);
        mutationalLoadChart.enabled(hasReliablePurityFit);
        mutationalLoadChart.scale(InlineBarChart.LOG10_SCALE);
        mutationalLoadChart.setTickMarks(new double[] { MutationalLoad.RANGE_MIN, 10, 100, MutationalLoad.RANGE_MAX }, noDecimalFormat);
        mutationalLoadChart.enableUndershoot(noDecimalFormat.format(0));
        mutationalLoadChart.enableOvershoot(">" + noDecimalFormat.format(mutationalLoadChart.max()));
        mutationalLoadChart.setIndicator(MutationalLoad.THRESHOLD,
                "Eligible for \nDRUP (" + noDecimalFormat.format(MutationalLoad.THRESHOLD) + ")");

        reportDocument.add(createCharacteristicDiv("Tumor mutational load",
                mutationalLoadString,
                "The tumor mutational load represents the total number of somatic missense variants across "
                        + "the whole genome of the tumor. Patients with a mutational load over 140 could be eligible for "
                        + "immunotherapy within the DRUP study.",
                mutationalLoadChart));

        double mutationalBurden = patientReport.tumorMutationalBurden();
        String mutationalBurdenString =
                hasReliablePurityFit ? singleDecimalFormat.format(mutationalBurden) + " variants per Mb" : DataUtil.NA_STRING;
        BarChart mutationalBurdenChart =
                new BarChart(mutationalBurden, MutationalBurden.RANGE_MIN, MutationalBurden.RANGE_MAX, "Low", "High", false);
        mutationalBurdenChart.enabled(hasReliablePurityFit);
        mutationalBurdenChart.scale(InlineBarChart.LOG10_SCALE);
        mutationalBurdenChart.setTickMarks(new double[] { MutationalBurden.RANGE_MIN, 10, MutationalBurden.RANGE_MAX },
                doubleDecimalFormat);
        mutationalBurdenChart.enableUndershoot(noDecimalFormat.format(0));
        mutationalBurdenChart.enableOvershoot(">" + singleDecimalFormat.format(mutationalBurdenChart.max()));
        reportDocument.add(createCharacteristicDiv("Tumor mutational burden",
                mutationalBurdenString,
                "The tumor mutational burden score represents the number of all somatic variants across the "
                        + "whole genome of the tumor per Mb.",
                mutationalBurdenChart));
    }

    @NotNull
    private Div createCharacteristicDiv(@NotNull String title, @NotNull String highlight, @NotNull String description,
            @NotNull BarChart chart) {
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

        return div;
    }
}
