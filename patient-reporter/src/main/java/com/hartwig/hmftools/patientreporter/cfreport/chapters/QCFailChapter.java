package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
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
        reportDocument.add(LineDivider.createLineDivider(getContentWidth()));

        // Content body
        LimsSampleType type = LimsSampleType.fromSampleId(failReport.sampleReport().sampleId());
        if (type == LimsSampleType.CORE) {
            reportDocument.add(createCoreContentBody());
        } else {
            reportDocument.add(createCPCTDRUPContentBody());
        }

        // End of report
        reportDocument.add(ReportSignature.createEndOfReportIndication().setMarginTop(20));
        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setMarginTop(20));

    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason, @NotNull SampleReport sampleReport) {

        if (sampleReport.recipient() == null) {
            throw new IllegalStateException("No recipient address present for sample " + sampleReport.sampleId());
        }

        final String title;
        final String reason;
        final String explanation;

        switch (failReason) {
            case LOW_DNA_YIELD: {
                title = "Notification tumor sample on hold for sequencing";
                reason = "Insufficient amount of DNA";
                explanation = "The amount of isolated DNA was <50 ng, which is insufficient for sequencing. ";
                break;
            }
            case LOW_TUMOR_PERCENTAGE: {
                title = "Notification of inadequate tumor sample";
                reason = "Not enough tumor cells detected by Pathology UMC Utrecht.";
                explanation = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            case POST_ANALYSIS_FAIL: {
                title = "Notification of inadequate tumor sample";
                reason = "Analysis has failed post DNA isolation";
                explanation = "This sample could not be processed to completion successfully.";
                break;
            }
            case SHALLOW_SEQ_LOW_PURITY: {
                title = "Notification of inadequate tumor sample";
                reason = "Not enough tumor DNA detected by molecular T % estimate.";
                explanation = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            default: {
                title = "TITLE";
                reason = "SUB_TITLE";
                explanation = "MESSAGE";
            }
        }

        // Initialize div
        Div div = new Div();
        div.setKeepTogether(true);

        // Add title
        div.add(new Paragraph(title)
                .addStyle(ReportResources.smallBodyHeadingStyle()));
        div.add(new Paragraph(reason).addStyle(ReportResources.dataHighlightStyle()));
        div.add(new Paragraph(explanation)
                .addStyle(ReportResources.bodyTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING));

        return div;

    }

    @NotNull
    private Div createCoreContentBody() {
        return createContentBody(new String[]{
                notSequencedText(),
                "When possible, please resubmit using the same DVO with project name " + failReport.sampleReport().projectName() + ".",
                "The HMF sample ID is " + failReport.sampleReport().sampleId() + " and the hospital patient ID is "
                        + failReport.sampleReport().hospitalPatientId(),
                "The project name of sample is " + failReport.sampleReport().projectName() + " and the submission ID is "
                        + failReport.sampleReport().submission(),
                "The internal tumor barcode is " + failReport.sampleReport().barcodeTumor() + " and the internal blood barcode is "
                        + failReport.sampleReport().barcodeReference(),
                "The tumor percentage estimated by Pathology UMC Utrecht is " + failReport.sampleReport().pathologyTumorPercentage(),
                shallowSeqText(),
                sampleArrivalDateText(),
                recipientText(),
                "The contact details are : " + failReport.sampleReport().contactNames() + " (" + failReport.sampleReport().contactEmails() + ")",
                accreditationText(),
                questionsText()

        });
    }

    @NotNull
    private Div createCPCTDRUPContentBody() {
        return createContentBody(new String[]{
                notSequencedText(),
                "When possible, please resubmit using the same " + failReport.study().studyName() + "-number. "
                        + "In case additional tumor material cannot be provided, please be notified that the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study.",
                "The HMF sample ID is " + failReport.sampleReport().sampleId(),
                "The internal tumor barcode is " + failReport.sampleReport().barcodeTumor() + " and the internal blood barcode is "
                        + failReport.sampleReport().barcodeReference(),
                "The tumor percentage estimated by Pathology UMC Utrecht is " + failReport.sampleReport().pathologyTumorPercentage(),
                shallowSeqText(),
                sampleArrivalDateText(),
                recipientText(),
                accreditationText(),
                questionsText()
        });
    }

    @NotNull
    private Div createContentBody(@NotNull String[] content) {
        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(getContentWidth());
        for (String line: content) {
            div.add(createContentParagraph(line));
        }
        return div;
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String text) {
        return new Paragraph(text)
                .addStyle(ReportResources.smallBodyTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static String questionsText() {
        return "For questions, please contact us via info@hartwigmedicalfoundation.nl";
    }

    @NotNull
    private String shallowSeqText() {
        return "The tumor percentage estimated by molecular tumor percentage is  " + failReport.sampleReport().purityShallowSeq();
    }

    @NotNull
    private String sampleArrivalDateText() {
        return "The biopsies evaluated for this sample have arrived on " + DataUtil.formatDate(failReport.sampleReport().tumorArrivalDate()) + " at "
                        + ReportResources.HARTWIG_ADDRESS;
    }

    @NotNull
    private static String notSequencedText() {
        return "The received biopsies for the tumor sample for this patient were inadequate to obtain a reliable sequencing "
                + "result. Therefore whole genome sequencing cannot be performed, "
                + "unless additional fresh tumor material can be provided for a new assessment.";
    }

    @NotNull
    private String recipientText() {
        return "This report is generated and verified by: " + failReport.user() + " and is addressed at " + failReport.sampleReport().recipient();
    }

    @NotNull
    private static String accreditationText() {
        return "The results on this report are based on tests that are performed under ISO/ICE-17025:2005 accreditation.";
    }

}
