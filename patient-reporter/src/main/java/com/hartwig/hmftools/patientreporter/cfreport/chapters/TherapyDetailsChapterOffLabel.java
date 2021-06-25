package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public class TherapyDetailsChapterOffLabel implements ReportChapter {

    @NotNull
    private final GenomicAnalysis genomicAnalysis;

    public TherapyDetailsChapterOffLabel(@NotNull final GenomicAnalysis genomicAnalysis) {
        this.genomicAnalysis = genomicAnalysis;
    }

    @NotNull
    @Override
    public String name() {
        return "Therapy details (Other tumor types)";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        Table chapterTable = new Table(1);

        chapterTable.addCell(new Cell().add(TherapyDetailsChapterFunctions.createEvidenceTable("Evidence on other tumor types",
                genomicAnalysis.offLabelEvidence())).setPadding(0).setBorder(Border.NO_BORDER));

        chapterTable.addFooterCell(new Cell().add(TherapyDetailsChapterFunctions.createChapterFootnoteOffLabel())
                .setPadding(0)
                .setBorder(Border.NO_BORDER));

        reportDocument.add(chapterTable);
    }
}
