package com.hartwig.hmftools.patientreporter.cfreport.components.tables;

public class EvidenceTable extends ReportTable {

    private static final float[] columnWidths = {18, 12, 61, 8, 15, 6};
    private static String[] headerNames = {"Drivers", "Match", "Treatments", "Level of evidence", "Response", "Source"};

    public EvidenceTable() {
        super(columnWidths, headerNames);
    }

    public void addRow(String driver, String match, String treatments, String evidenceLevel, String response, String source) {
        addCell(driver);
        addCell(match);
        addCell(treatments);
        addCell(evidenceLevel);
        addCell(response);
        addCell(source);
    }


}
