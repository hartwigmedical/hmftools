package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.itextpdf.layout.Document;

public class DetailsAndDisclaimerChapter extends ReportChapter {
    @Override
    public String getName() {
        return "Sample details & disclaimer";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ClosingChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

    }
}
