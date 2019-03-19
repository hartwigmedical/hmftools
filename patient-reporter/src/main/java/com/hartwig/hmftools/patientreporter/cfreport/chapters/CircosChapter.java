package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

public class CircosChapter extends ReportChapter {
    @Override
    public String getName() {
        return "CIRCOS plot";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

    }
}
