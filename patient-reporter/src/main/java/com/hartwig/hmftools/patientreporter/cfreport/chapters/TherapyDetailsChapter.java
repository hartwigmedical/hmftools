package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.SectionTitle;
import com.hartwig.hmftools.patientreporter.cfreport.components.tables.ClinicalTrialsTable;
import com.hartwig.hmftools.patientreporter.cfreport.components.tables.EvidenceTable;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
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

        report.add(SectionTitle.getSectionTitle("Tumor type specific evidence"));

        EvidenceTable table = new EvidenceTable();
        for (int i = 0; i < 10; i++) {
            table.addRow("BRAF p.Val600Glu", "Specific", "Binimetinib + Encorafenib", "A", "Responsive", "OncoKB");
        }
        report.add(table.getTable());

    }

    private static void addClinicalTrials(@NotNull Document report) {

        report.add(SectionTitle.getSectionTitle("Clinical trials (NL)"));

        ClinicalTrialsTable table = new ClinicalTrialsTable();
        for (int i = 0; i < 10; i++) {
            table.addRow("BRAF p.Val600Glu", "Specific", "Binimetinib + Encorafenib", "NL57739.031.16", "IClusion");
        }
        report.add(table.getTable());

    }

    private static void addOtherTumorTypeEvidence(@NotNull Document report) {

        report.add(SectionTitle.getSectionTitle("Evidence on other tumor types"));

        EvidenceTable table = new EvidenceTable();
        for (int i = 0; i < 50; i++) {
            table.addRow("BRAF p.Val600Glu", "Specific", "Binimetinib + Encorafenib", "A", "Responsive", "OncoKB");
        }
        report.add(table.getTable());

    }


}
