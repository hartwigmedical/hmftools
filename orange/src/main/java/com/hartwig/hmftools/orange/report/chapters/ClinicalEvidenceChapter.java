package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class ClinicalEvidenceChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public ClinicalEvidenceChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Clinical Evidence";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Clinical Evidence"));
        document.add(new Paragraph("TODO: All on-label evidence, grouped by treatment").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: All off-label evidence, grouped by treatment").addStyle(ReportResources.tableContentStyle()));
    }
}
