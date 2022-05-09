package com.hartwig.hmftools.patientreporter.cfreport.chapters.panel;

import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReason;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class PanelQCFailChapter implements ReportChapter {

    private static final String TITLE_REPORT = "Oncopanel Result Report Failed";

    @NotNull
    private final PanelFailReport report;

    public PanelQCFailChapter(@NotNull PanelFailReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return report.isCorrectedReport() ? TITLE_REPORT + " (Corrected)" : TITLE_REPORT;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return TITLE_REPORT;
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        reportDocument.add(TumorLocationAndTypeTable.createBiopsyLocationAndTumorLocation(report.sampleReport()
                .primaryTumorLocationString(), report.sampleReport().biopsyLocationString(), contentWidth()));
        reportDocument.add(new Paragraph());
        reportDocument.add(TumorLocationAndTypeTable.createTumorType(report.sampleReport().primaryTumorTypeString(), contentWidth()));

        reportDocument.add(new Paragraph("The information regarding 'primary tumor location', 'primary tumor type' and 'biopsy location'"
                + " is based on information received from the originating hospital.").addStyle(ReportResources.subTextSmallStyle()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));
        reportDocument.add(createFailReasonDiv(report.panelFailReason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));
    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull PanelFailReason failReason) {
        String reason = DataUtil.NA_STRING;
        String explanation = DataUtil.NA_STRING;
        String explanationDetail = DataUtil.NA_STRING;

        switch (failReason) {
            case PANEL_FAILURE: {
                reason = "Insufficient quality of received biomaterial(s)";
                explanation = "The received biomaterial(s) did not meet the requirements that are needed for \n"
                        + "high quality Next Generation Sequencing.";
                explanationDetail =
                        "Sequencing could not be performed due to insufficient DNA.";
                break;
            }
        }

        Div div = new Div();
        div.setKeepTogether(true);

        div.add(new Paragraph(reason).addStyle(ReportResources.dataHighlightStyle()));
        div.add(new Paragraph(explanation).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING));
        div.add(new Paragraph(explanationDetail).addStyle(ReportResources.subTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING));
        return div;
    }
}