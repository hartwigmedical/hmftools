package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

public class ExplanationChapter extends ReportChapter {
    @Override
    public String getName() {
        return "Report explanation";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

    }
}
