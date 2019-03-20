package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.components.SectionTitle;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Table;
import org.jetbrains.annotations.NotNull;

public class TherapyDetailsChapter extends ReportChapter {

    @Override
    public String getName() {
        return "Therapy details";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {
        addTumorTypeEvidence(report);
        addClinicalTrials(report);
        addOtherTumorTypeEvidence(report);
    }

    private static void addTumorTypeEvidence(@NotNull Document report) {

        report.add(SectionTitle.createSectionTitle("Tumor type specific evidence"));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 61, 8, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "Level of evidence", "Response", "Source"});

        for (int i = 0; i < 10; i++) {
            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell("Specific"));
            table.addCell(TableHelper.getContentCell("Binimetinib + Encorafenib"));
            table.addCell(TableHelper.getContentCell("A"));
            table.addCell(TableHelper.getContentCell("Responsive"));
            table.addCell(TableHelper.getContentCell("OncoKB"));
        }

        report.add(table);

    }

    private static void addClinicalTrials(@NotNull Document report) {

        report.add(SectionTitle.createSectionTitle("Clinical trials (NL)"));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 69, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "CCMO", "Source"});

        for (int i = 0; i < 10; i++) {
            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell("Specific"));
            table.addCell(TableHelper.getContentCell("Binimetinib + Encorafenib"));
            table.addCell(TableHelper.getContentCell("NL57739.031.16"));
            table.addCell(TableHelper.getContentCell("IClusion"));
        }

        report.add(table);

    }

    private static void addOtherTumorTypeEvidence(@NotNull Document report) {

        report.add(SectionTitle.createSectionTitle("Evidence on other tumor types"));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 61, 8, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "Level of evidence", "Response", "Source"});

        for (int i = 0; i < 10; i++) {
            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell("Specific"));
            table.addCell(TableHelper.getContentCell("Binimetinib + Encorafenib"));
            table.addCell(TableHelper.getContentCell("A"));
            table.addCell(TableHelper.getContentCell("Responsive"));
            table.addCell(TableHelper.getContentCell("OncoKB"));
        }

        report.add(table);

    }


}
