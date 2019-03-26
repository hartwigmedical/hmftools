package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.hartwig.hmftools.patientreporter.cfreport.components.Icon;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.*;
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

        Table chapterTable = new Table(1);

        chapterTable.addCell(new Cell()
                .add(createTumorTypeEvidenceTable())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addCell(new Cell()
                .add(createClinicalTrialsTable())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        chapterTable.addCell(new Cell()
                .add(createOtherTumorTypeEvidenceTable())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        // Add legend/footnote that will appear on each page of the chapter
        chapterTable.addFooterCell(new Cell()
                .add(createChapterFootnote())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        report.add(chapterTable);

    }

    @NotNull
    private static Table createTumorTypeEvidenceTable() {

        // Temporary content
        // @TODO remove
        final String[] levels = {"A", "B", "C", "D"};
        final String[][] treatmentCombinations = {
                {"Binimetinib", "Encorafenib" },
                {"RO4987655"},
                {"Dabrafenib"},
                {"Dabrafenib", "Trametinib"},
                {"Vemurafenib"},
                {"Cobimetinib", "Vemurafenib", "Encorafenib"},
                {"Encorafenib"}
        };

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {18, 12, 61, 8, 15, 6});
        for (int i = 0; i < 10; i++) {

            final String level = levels[(int) (Math.random() * (float) levels.length)];
            String[] treatments = treatmentCombinations[(int) (Math.random() * (float) treatmentCombinations.length)];

            contentTable.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentParagraph(treatments)).setVerticalAlignment(VerticalAlignment.TOP));
            contentTable.addCell(TableHelper.getContentCell(new Paragraph(Icon.createLevelIcon(level))));
            contentTable.addCell(TableHelper.getContentCell("Responsive"));
            contentTable.addCell(TableHelper.getContentCell("OncoKB"));

        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable("Tumor type specific evidence", contentTable, new Cell[]  {
                TableHelper.getHeaderCell("Drivers"),
                TableHelper.getHeaderCell("Match"),
                TableHelper.getHeaderCell("Treatments"),
                TableHelper.getHeaderCell("Level of evidence"),
                TableHelper.getHeaderCell("Response"),
                TableHelper.getHeaderCell("Source")
        });

    }

    @NotNull
    private static Table createClinicalTrialsTable() {

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {18, 12, 69, 15, 6});
        for (int i = 0; i < 20; i++) {
            contentTable.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentParagraph(new String[] {"Binimetinib", "Encorafenib"})));
            contentTable.addCell(TableHelper.getContentCell("NL57739.031.16"));
            contentTable.addCell(TableHelper.getContentCell("IClusion"));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable("Clinical trials (NL)", contentTable, new Cell[] {
                TableHelper.getHeaderCell("Drivers"),
                TableHelper.getHeaderCell("Match"),
                TableHelper.getHeaderCell("Treatments"),
                TableHelper.getHeaderCell("CCMO"),
                TableHelper.getHeaderCell("Source")
        });

    }

    @NotNull
    private static Table createOtherTumorTypeEvidenceTable() {

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {18, 12, 61, 8, 15, 6});
        for (int i = 0; i < 10; i++) {
            contentTable.addCell(TableHelper.getContentCell("BRAF p.Val600Glu"));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentMatchParagraph(Math.random() > 0.5)));
            contentTable.addCell(TableHelper.getContentCell(createTreatmentParagraph(new String[] {"Vemurafenib"})));
            contentTable.addCell(TableHelper.getContentCell(new Paragraph(Icon.createLevelIcon("A"))));
            contentTable.addCell(TableHelper.getContentCell("Responsive"));
            contentTable.addCell(TableHelper.getContentCell("OncoKB"));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable("Evidence on other tumor types", contentTable, new Cell[] {
                TableHelper.getHeaderCell("Drivers"),
                TableHelper.getHeaderCell("Match"),
                TableHelper.getHeaderCell("Treatments"),
                TableHelper.getHeaderCell("Level of evidence"),
                TableHelper.getHeaderCell("Response"),
                TableHelper.getHeaderCell("Source")
        });

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

    @NotNull
    private final static Paragraph createChapterFootnote() {
        return new Paragraph()
                .setKeepTogether(true)
                .add("The Cancer Genome Interpreter (CGI), OncoKB and CiViC knowledge bases are used to " +
                    "annotate variants of all types with clinical evidence. Only treatment associated evidence with a high " +
                    "level of evidence ( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_A))
                .add(" validated association; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_B))
                .add(" strong clinical evidence) are reported here. Potential evidence items with a lower level of evidence ( ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_C))
                .add(" case study, limited clinical evidence; ")
                .add(Icon.createIcon(Icon.IconType.LEVEL_D))
                .add(" pre-clinical) are not reported.")
                    .addStyle(ReportResources.subTextStyle());
    }

}
