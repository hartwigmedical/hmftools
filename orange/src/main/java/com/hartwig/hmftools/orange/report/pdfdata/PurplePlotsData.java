package com.hartwig.hmftools.orange.report.pdfdata;

public class PurplePlotsData
{
    public final boolean hasPurpleFail;
    public final String purpleInputCircosPlotPath;
    public final String purpleCopyNumberPlotPath;
    public final String purpleClonalityPlotPath;
    public final String purplePurityRangePlotPath;
    public final String purpleMinorAlleleMapPlotPath;
    public final String purpleRainfallPlotPath;

    public PurplePlotsData(
            final boolean hasPurpleFail,
            final String purpleInputCircosPlotPath,
            final String purpleCopyNumberPlotPath,
            final String purpleClonalityPlotPath,
            final String purplePurityRangePlotPath,
            final String purpleMinorAlleleMapPlotPath,
            final String purpleRainfallPlotPath)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.purpleInputCircosPlotPath = purpleInputCircosPlotPath;
        this.purpleCopyNumberPlotPath = purpleCopyNumberPlotPath;
        this.purpleClonalityPlotPath = purpleClonalityPlotPath;
        this.purplePurityRangePlotPath = purplePurityRangePlotPath;
        this.purpleMinorAlleleMapPlotPath = purpleMinorAlleleMapPlotPath;
        this.purpleRainfallPlotPath = purpleRainfallPlotPath;
    }
}
