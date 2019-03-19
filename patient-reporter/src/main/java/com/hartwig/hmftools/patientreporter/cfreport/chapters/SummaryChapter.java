package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

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

    }

}
