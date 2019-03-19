package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

public class SummaryChapter extends ReportChapter {

    @Override
    public String getName() {
        return "Summary";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.SummaryChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {
        report.add(new Paragraph("Summary content"));

    }

}
