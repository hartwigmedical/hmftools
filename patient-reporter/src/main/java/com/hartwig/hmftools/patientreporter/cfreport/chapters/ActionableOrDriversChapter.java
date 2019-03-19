package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

public class ActionableOrDriversChapter extends ReportChapter {

    @Override
    public String getName() {
        return "Actionable or drivers";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

    }

}
