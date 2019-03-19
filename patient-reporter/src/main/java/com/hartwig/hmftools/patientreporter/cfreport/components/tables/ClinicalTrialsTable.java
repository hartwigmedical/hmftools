package com.hartwig.hmftools.patientreporter.cfreport.components.tables;

public class ClinicalTrialsTable extends ReportTable {

    private static final float[] columnWidths = {18, 12, 69, 15, 6};
    private static String[] headerNames = {"Drivers", "Match", "Treatments", "CCMO", "Source"};

    public ClinicalTrialsTable() {
        super(columnWidths, headerNames);
    }

    public void addRow(String driver, String match, String treatments, String ccmo, String source) {
        addCell(driver);
        addCell(match);
        addCell(treatments);
        addCell(ccmo);
        addCell(source);
    }
}
