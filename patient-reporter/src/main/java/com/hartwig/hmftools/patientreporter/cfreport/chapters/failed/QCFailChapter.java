package com.hartwig.hmftools.patientreporter.cfreport.chapters.failed;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class QCFailChapter implements ReportChapter {

    private static final String TITLE_REPORT = "Failed DNA Analysis Report";

    @NotNull
    private final QCFailReport failReport;

    public QCFailChapter(@NotNull QCFailReport failReport) {
        this.failReport = failReport;
    }

    @NotNull
    @Override
    public String name() {
        return failReport.isCorrectedReport() ? TITLE_REPORT + " (Corrected)" : TITLE_REPORT;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
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
        reportDocument.add(TumorLocationAndTypeTable.createBiopsyLocationAndTumorLocation(failReport.sampleReport()
                .primaryTumorLocationString(), failReport.sampleReport().biopsyLocationString(), contentWidth()));
        reportDocument.add(new Paragraph());
        reportDocument.add(TumorLocationAndTypeTable.createTumorType(failReport.sampleReport().primaryTumorTypeString(),
                contentWidth()));

        reportDocument.add(new Paragraph("The information regarding 'primary tumor location', 'primary tumor type' and 'biopsy location'"
                + " is based on information received from the originating hospital.").addStyle(ReportResources.subTextSmallStyle()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason) {
        String reason = DataUtil.NA_STRING;
        String explanation = DataUtil.NA_STRING;
        String explanationDetail = DataUtil.NA_STRING;
        boolean reportPeachReport = true;

        switch (failReason.type()) {
            case LOW_QUALITY_BIOPSY: {
                reason = "Insufficient quality of received biomaterial(s)";
                explanation = "The received biomaterial(s) did not meet the requirements that are needed for \n"
                        + "high quality Whole Genome Sequencing.";
                break;
            }
            case TECHNICAL_FAILURE: {
                reason = "Technical failure";
                explanation = "Whole Genome Sequencing could not be successfully performed on the received biomaterial(s) \n"
                        + "due to technical problems.";
            }
        }

        switch (failReason) {
            case INSUFFICIENT_DNA: {
                explanationDetail =
                        "Sequencing could not be performed due to insufficient DNA.";
                break;
            }
            case INSUFFICIENT_TCP_DEEP_WGS:
            case INSUFFICIENT_TCP_SHALLOW_WGS: {
                explanationDetail = "The tumor percentage based on molecular estimation was below the minimal of 20% tumor cells \n"
                        + "and could not be further analyzed.";
                break;
            }
            case SUFFICIENT_TCP_QC_FAILURE: {
                explanationDetail = "The tumor percentage based on molecular estimation was above the minimal of 20% tumor cells \n"
                        + "but could not be further analyzed due to insufficient quality.";
                break;
            }
            case TECHNICAL_FAILURE: {
                explanationDetail = Strings.EMPTY;
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
