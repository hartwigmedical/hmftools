package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

public class TumorCharacteristicsChapter extends ReportChapter {

    @Override
    public String getName() {
        return "Tumor characteristics";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

    }
}
