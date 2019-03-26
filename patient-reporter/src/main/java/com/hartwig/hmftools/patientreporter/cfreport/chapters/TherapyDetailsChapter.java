package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

import java.util.StringJoiner;

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

        report.add(new Paragraph("Tumor type specific evidence")
                .addStyle(ReportResources.sectionTitleStyle()));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 61, 8, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "Level of evidence", "Response", "Source"});


        String[] levels = {"A", "B", "C", "D"};
        String[][] treatmentCombinations = {
                {"Binimetinib", "Encorafenib" },
                {"RO4987655"},
                {"Dabrafenib"},
                {"Dabrafenib", "Trametinib"},
                {"Vemurafenib"},
                {"Cobimetinib", "Vemurafenib", "Encorafenib"},
                {"Encorafenib"}
        };
        for (int i = 0; i < 10; i++) {

            final String level = levels[(int) (Math.random() * (float) levels.length)];
            String[] treatments = treatmentCombinations[(int) (Math.random() * (float) treatmentCombinations.length)];

            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            table.addCell(TableHelper.getContentCell(createTreatmentParagraph(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            table.addCell(TableHelper.getContentCell(new Paragraph(Icon.createLevelIcon(level))));
            table.addCell(TableHelper.getContentCell("Responsive"));
            table.addCell(TableHelper.getContentCell("OncoKB"));
        }

        report.add(table);

    }

    private static void addClinicalTrials(@NotNull Document report) {

        report.add(new Paragraph("Clinical trials (NL)")
                .addStyle(ReportResources.sectionTitleStyle()));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 69, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "CCMO", "Source"});

        for (int i = 0; i < 10; i++) {
            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            table.addCell(TableHelper.getContentCell(createTreatmentParagraph(new String[] {"Binimetinib", "Encorafenib"})));
            table.addCell(TableHelper.getContentCell("NL57739.031.16"));
            table.addCell(TableHelper.getContentCell("IClusion"));
        }

        report.add(table);

    }

    private static void addOtherTumorTypeEvidence(@NotNull Document report) {

        report.add(new Paragraph("Evidence on other tumor types")
                .addStyle(ReportResources.sectionTitleStyle()));

        Table table = TableHelper.getReportTable(
                new float[] {18, 12, 61, 8, 15, 6},
                new String[] {"Drivers", "Match", "Treatments", "Level of evidence", "Response", "Source"});

        for (int i = 0; i < 10; i++) {
            table.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            table.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            table.addCell(TableHelper.getContentCell(createTreatmentParagraph(new String[] {"Vemurafenib"})));
            table.addCell(TableHelper.getContentCell(new Paragraph(Icon.createLevelIcon("A"))));
            table.addCell(TableHelper.getContentCell("Responsive"));
            table.addCell(TableHelper.getContentCell("OncoKB"));
        }

        report.add(table);

    }

    @NotNull
    private final static Paragraph createTreatmentMatchParagraph(boolean matchIsSpecific) {
        return new Paragraph()
                .add(Icon.createIcon(matchIsSpecific ? Icon.IconType.MATCH_SPECIFIC : Icon.IconType.MATCH_BROAD))
                .add(new Text(" " + (matchIsSpecific ? "Specific" : "Broad"))
                        .addStyle(ReportResources.tableContentStyle()));
    }

    @NotNull
    private final static Paragraph createTreatmentParagraph(@NotNull String[] treatments) {

        Paragraph p = new Paragraph();
        StringJoiner joiner = new StringJoiner(" + ");
        for (String treatmentName: treatments) {
            treatmentName = treatmentName.trim();

            // Add treatment icon
            p.add(Icon.createTreatmentIcon(treatmentName));

            // Add treatment name to joiner
            joiner.add(treatmentName);

        }

        // Add treatment names
        p.add(new Text(" " + joiner.toString())
                        .addStyle(ReportResources.tableContentStyle()));

        return p;
    }

}
