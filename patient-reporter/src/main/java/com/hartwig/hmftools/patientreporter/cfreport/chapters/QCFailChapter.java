package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class QCFailChapter implements ReportChapter {

    @NotNull
    private final QCFailReport failReport;

    public QCFailChapter(@NotNull QCFailReport failReport) {
        this.failReport = failReport;
    }

    @NotNull
    @Override
    public String name() {
        return failReport.reason().title();
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
        if (failReport.sampleReport().addressee() == null) {
            throw new IllegalStateException("No recipient address present for sample " + failReport.sampleReport().sampleId());
        }

        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(failReport.sampleReport().primaryTumorLocationString(),
                failReport.sampleReport().cancerSubTypeString(),
                contentWidth()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()).setMarginTop(10));

        LimsSampleType type = LimsSampleType.fromSampleId(failReport.sampleReport().sampleId());

        switch (type) {
            case CORE:
                reportDocument.add(createCoreContentBody());
                break;
            case WIDE:
                reportDocument.add(createWIDEContentBody());
                break;
            default:
                reportDocument.add(createCPCTDRUPContentBody());
        }

        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setMarginTop(20));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason) {
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
                reason = "Not enough tumor DNA detected based on molecular estimate.";
                explanation = "For sequencing we require a minimum of 20% tumor cells.";
                break;
            }
            case INSUFFICIENT_TISSUE: {
                title = "Notification tumor sample on hold for sequencing";
                reason = "Insufficient tissue";
                explanation = "The amount of isolated DNA was <50 ng, which is insufficient for sequencing. ";
                break;
            }
            default: {
                title = "TITLE";
                reason = "SUB_TITLE";
                explanation = "MESSAGE";
            }
        }

        Div div = new Div();
        div.setKeepTogether(true);

        div.add(new Paragraph(title.toUpperCase()).addStyle(ReportResources.subTextStyle()));
        div.add(new Paragraph(reason).addStyle(ReportResources.dataHighlightStyle()));
        div.add(new Paragraph(explanation).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING));

        return div;
    }

    @NotNull
    private Table createWIDEContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createWideContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createWideContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createWideContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph(
                "If available new tumor material can be provided for a new assessment, please resubmit using the same " + failReport.study()
                        .studyName() + "-number. " + "If additional material cannot be provided the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study."));
        divColumn1.add(createContentParagraphTwice("The HMF sample ID is ",
                failReport.sampleReport().sampleId(),
                " and the tissue ID of pathology is: ",
                failReport.sampleReport().hospitalPathologySampleId()));
        divColumn1.add(obtainedResults());
        divColumn1.add(tumorSampleDataText());
        divColumn1.add(bloodSampleDataText());
        divColumn1.add(relatedSamples());
        divColumn1.add((shallowSeqText()));
        return divColumn1;
    }

    @NotNull
    private Div createWideContentBodyColumn2() {
        Div divColumn2 = new Div();

        divColumn2.add((evaluatedAddress()));
        divColumn2.add((recipientText()));
        divColumn2.add(versionPatientReport());
        divColumn2.add((accreditationText()));
        divColumn2.add(disclaimerTumorLocation());
        divColumn2.add(sentivityResults());
        divColumn2.add((questionsText()));
        return divColumn2;
    }

    @NotNull
    private Table createCoreContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCoreContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createCoreContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createCoreContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph("If available new tumor material can be provided for a new assessment, "
                + "please resubmit using the same DVO with project name " + failReport.sampleReport().projectName() + "."));
        divColumn1.add(createContentParagraphTwice("The HMF sample ID is ",
                failReport.sampleReport().sampleId(),
                " and the hospital patient ID is ",
                failReport.sampleReport().hospitalPatientId()));
        divColumn1.add(createContentParagraphTwice("The project name of sample is ",
                failReport.sampleReport().projectName(),
                " and the submission ID is ",
                failReport.sampleReport().submissionId()));
        divColumn1.add(obtainedResults());
        divColumn1.add(tumorSampleDataText());
        divColumn1.add(bloodSampleDataText());
        divColumn1.add(relatedSamples());
        divColumn1.add(createContentParagraph("The tumor percentage estimated by Pathology UMC Utrecht is ",
                failReport.sampleReport().pathologyTumorPercentage()));
        divColumn1.add((shallowSeqText()));
        return divColumn1;
    }

    @NotNull
    private Div createCoreContentBodyColumn2() {
        Div divColumn2 = new Div();
        divColumn2.add((evaluatedAddress()));
        divColumn2.add((recipientText()));
        divColumn2.add(createContentParagraphRequest(failReport.sampleReport()));
        divColumn2.add(versionPatientReport());
        divColumn2.add(accreditationText());
        divColumn2.add(disclaimerTumorLocation());
        divColumn2.add(sentivityResults());
        divColumn2.add(questionsText());
        return divColumn2;
    }

    @NotNull
    private static Paragraph createContentParagraphRequest(@NotNull SampleReport sampleReport) {
        return createContentParagraph("The requester is : ").add(new Text(sampleReport.requesterName()).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add("(")
                .add(new Text(sampleReport.requesterEmail()).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(")")
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private Table createCPCTDRUPContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCPCTDRUPContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createCPCTDRUPContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createCPCTDRUPContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph(
                "If available new tumor material can be provided for a new assessment, please resubmit using the same " + failReport.study()
                        .studyName() + "-number. " + "If additional material cannot be provided the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study."));
        divColumn1.add(createContentParagraph("The HMF sample ID is ", failReport.sampleReport().sampleId()));
        divColumn1.add(obtainedResults());
        divColumn1.add(tumorSampleDataText());
        divColumn1.add(bloodSampleDataText());
        divColumn1.add(createContentParagraph("The tumor percentage estimated by Pathology UMC Utrecht is ",
                failReport.sampleReport().pathologyTumorPercentage()));
        divColumn1.add(relatedSamples());
        divColumn1.add(shallowSeqText());
        return divColumn1;
    }

    @NotNull
    private Div createCPCTDRUPContentBodyColumn2() {
        Div divColumn2 = new Div();
        divColumn2.add(evaluatedAddress());
        divColumn2.add(recipientText());
        divColumn2.add(versionPatientReport());
        divColumn2.add(accreditationText());
        divColumn2.add(disclaimerTumorLocation());
        divColumn2.add(sentivityResults());
        divColumn2.add(questionsText());
        return divColumn2;
    }

    @NotNull
    private static Paragraph relatedSamples() {
        return createContentParagraph("The results stated in these report are based on the tested tumor and blood sample.");
    }

    @NotNull
    private Paragraph obtainedResults() {
        String earliestArrivalDate = failReport.sampleReport().earliestArrivalDate();
        return createContentParagraphTwice("The results in this report have been obtained between ",
                earliestArrivalDate != null ? earliestArrivalDate : DataUtil.NA_STRING,
                " and ",
                ReportResources.REPORT_DATE);
    }

    @NotNull
    private Paragraph sentivityResults() {
        return createContentParagraph("Based on a tumor purity of at least 30%, the test has a sensitivity of >95% for detection "
                + "of somatic variants and >95% for detection of translocations and gene copy number changes. For samples with a "
                + "purity above 20%, the test has a sensitivity of >90%.");
    }

    @NotNull
    private Paragraph disclaimerTumorLocation() {
        return createContentParagraph("The ‘primary tumor location’ and ‘cancer subtype’ are received from the requesting hospital and "
                + "have influence on the clinical evidence/study matching. No check is performed to verify the received information.");
    }

    @NotNull
    private static Paragraph versionPatientReport() {
        return createContentParagraph("This report is based by patient reporter ", ReportResources.VERSION_REPORT);
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph questionsText() {
        return createContentParagraph("For questions, please contact us via ", "info@hartwigmedicalfoundation.nl");
    }

    @NotNull
    private Paragraph shallowSeqText() {
        return createContentParagraph("The tumor percentage based on molecular estimation is ",
                failReport.sampleReport().purityShallowSeq());
    }

    @NotNull
    private Paragraph tumorSampleDataText() {
        return createContentParagraphTwice("This experiment is performed on the tumor sample which arrived on ",
                DataUtil.formatDate(failReport.sampleReport().tumorArrivalDate()),
                " with internal tumor barcode ",
                failReport.sampleReport().tumorBarcode());
    }

    @NotNull
    private Paragraph bloodSampleDataText() {
        return createContentParagraphTwice("This experiment is performed on the blood sample which arrived on ",
                DataUtil.formatDate(failReport.sampleReport().refArrivalDate()),
                " with internal blood barcode ",
                failReport.sampleReport().refBarcode());
    }

    @NotNull
    private Paragraph evaluatedAddress() {
        return createContentParagraph("The biopsies evaluated at ", ReportResources.HARTWIG_ADDRESS);
    }

    @NotNull
    private static Paragraph notSequencedText() {
        return createContentParagraph("The received tumor biopsies were inadequate for whole genome sequencing. ");
    }

    @NotNull
    private Paragraph recipientText() {
        String addressee = failReport.sampleReport().addressee();
        assert addressee != null; // Has been checked prior to calling this function.
        return createContentParagraph("This report is generated and verified by: " + failReport.user() + " and is addressed to ",
                addressee);
    }

    @NotNull
    private static Paragraph accreditationText() {
        return createContentParagraph(
                "The results on this report are based on tests that are performed under ISO/ICE-17025:2005 accreditation.");
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String regularPart, @NotNull String boldPart) {
        return createContentParagraph(regularPart).add(new Text(boldPart).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph createContentParagraphTwice(@NotNull String regularPart, @NotNull String boldPart,
            @NotNull String regularPart2, @NotNull String boldPart2) {
        return createContentParagraph(regularPart).add(new Text(boldPart).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(regularPart2)
                .add(new Text(boldPart2).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }
}
