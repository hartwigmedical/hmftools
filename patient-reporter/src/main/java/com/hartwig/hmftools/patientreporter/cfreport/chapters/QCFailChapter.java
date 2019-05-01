package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public class QCFailChapter implements ReportChapter {

    private final QCFailReport failReport;

    public QCFailChapter(QCFailReport failReport) {
        this.failReport = failReport;
    }

    @Override
    public String getName() {
        return failReport.reason().title();
    }

    public boolean isFullWidth() {
        return false;
    }

    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public void render(@NotNull Document reportDocument) throws IOException {

        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(
                failReport.sampleReport().primaryTumorLocationString(),
                failReport.sampleReport().cancerSubTypeString(),
                getContentWidth()));
        reportDocument.add(LineDivider.createLineDivider(getContentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason(), failReport.sampleReport()));

        // Content body
        LimsSampleType type = LimsSampleType.fromSampleId(failReport.sampleReport().sampleId());
        if (type == LimsSampleType.CORE) {
            reportDocument.add(createCoreContentBody());
        } else {
            reportDocument.add(createCPCTDRUPContentBody());
        }

        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setPaddingTop(80));

    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason, @NotNull SampleReport sampleReport) {

        if (sampleReport.recipient() == null) {
            throw new IllegalStateException("No recipient address present for sample " + sampleReport.sampleId());
        }

        final String title;
        final String subTitle;
        final String message;

        switch (failReason) {
            case LOW_DNA_YIELD: {
                title = "Notification tumor sample on hold for sequencing";
                subTitle = "Insufficient amount of DNA";
                message = "The amount of isolated DNA was <50 ng, which is insufficient for sequencing. ";
                break;
            }
            case LOW_TUMOR_PERCENTAGE: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Not enough tumor cells detected by Pathology UMC Utrecht.";
                message = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            case POST_ANALYSIS_FAIL: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Analysis has failed post DNA isolation";
                message = "This sample could not be processed to completion successfully.";
                break;
            }
            case SHALLOW_SEQ_LOW_PURITY: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Not enough tumor DNA detected by molecular T % estimate.";
                message = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            default: {
                title = "TITLE";
                subTitle = "SUB_TITLE";
                message = "MESSAGE";
            }
        }

        // Initialize div
        Div div = new Div();
        div.setKeepTogether(true);

        // Add title
        div.add(new Paragraph(title)
                .addStyle(ReportResources.sectionTitleStyle()));
        div.add(DataLabel.createDataLabel(subTitle));
        div.add(new Paragraph(message)
                .addStyle(ReportResources.bodyTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING));

        return div;

    }


    private static Div createCoreContentBody() {
        return new Div();
    }

    private static Div createCPCTDRUPContentBody() {
        return new Div();
    }

}
