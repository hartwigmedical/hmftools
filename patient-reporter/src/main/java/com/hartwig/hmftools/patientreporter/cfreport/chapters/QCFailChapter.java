package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
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
        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(failReport.sampleReport().primaryTumorLocationString(),
                failReport.sampleReport().cancerSubTypeString(),
                contentWidth()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason(), failReport.sampleReport()));
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

        reportDocument.add(ReportSignature.createEndOfReportIndication().setMarginTop(20));
        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setMarginTop(20));
    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason, @NotNull SampleReport sampleReport) {
        if (sampleReport.addressee() == null) {
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
                reason = "Not enough tumor DNA detected based on molecular estimate.";
                explanation = "For sequencing we require a minimum of 30% tumor cells.";
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
        table.addCell(TableUtil.createLayoutCell()); // Spacer
        table.addCell(TableUtil.createLayoutCell().add(createWideContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createWideContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph(
                "Please resubmit using the same " + failReport.study().studyName() + "-number. "
                        + "If additional material cannot be provided the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study."));
        divColumn1.add(createContentParagraphTwice(
                "The HMF sample ID is " , failReport.sampleReport().sampleId() , " and the tissue ID of pathology is: "
                        , failReport.sampleReport().hospitalPathologySampleId()));
        divColumn1.add(createContentParagraphTwice(
                "The internal tumor barcode is " , failReport.sampleReport().tumorBarcode() , " and the internal blood barcode is "
                        , failReport.sampleReport().refBarcode()));
        return divColumn1;
    }

    @NotNull
    private Div createWideContentBodyColumn2() {
        Div divColumn2 = new Div();

        divColumn2.add((shallowSeqText()));
        divColumn2.add((sampleArrivalDateText()));
        divColumn2.add((recipientText()));
        divColumn2.add((accreditationText()));
        divColumn2.add((questionsText()));
        return divColumn2;
    }

    @NotNull
    private Table createCoreContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCoreContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell()); // Spacer
        table.addCell(TableUtil.createLayoutCell().add(createCoreContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createCoreContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph(
                "Please resubmit using the same DVO with project name " + failReport.sampleReport().projectName() + "."));
        divColumn1.add(createContentParagraphTwice(
                "The HMF sample ID is " , failReport.sampleReport().sampleId() , " and the hospital patient ID is "
                        , failReport.sampleReport().hospitalPatientId()));
        divColumn1.add(createContentParagraphTwice(
                "The project name of sample is " , failReport.sampleReport().projectName() , " and the submission ID is "
                        , failReport.sampleReport().submissionId()));
        divColumn1.add(createContentParagraphTwice(
                "The internal tumor barcode is " , failReport.sampleReport().tumorBarcode() , " and the internal blood barcode is "
                        , failReport.sampleReport().refBarcode()));
        return divColumn1;
    }

    @NotNull
    private Div createCoreContentBodyColumn2() {
        Div divColumn2 = new Div();

        divColumn2.add(createContentParagraph(
                "The tumor percentage estimated by Pathology UMC Utrecht is " , failReport.sampleReport().pathologyTumorPercentage()));
        divColumn2.add((shallowSeqText()));
        divColumn2.add((sampleArrivalDateText()));
        divColumn2.add((recipientText()));
        divColumn2.add(createContentParagraphRequest(failReport.sampleReport()));
        divColumn2.add(accreditationText());
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
        table.addCell(TableUtil.createLayoutCell()); // Spacer
        table.addCell(TableUtil.createLayoutCell().add(createCPCTDRUPContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createCPCTDRUPContentBodyColumn1() {
        Div divColumn1 = new Div();
        divColumn1.add(notSequencedText());
        divColumn1.add(createContentParagraph(
                "Please resubmit using the same " + failReport.study().studyName() + "-number. "
                        + "If additional material cannot be provided the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study."));
        divColumn1.add(createContentParagraph("The HMF sample ID is " , failReport.sampleReport().sampleId()));
        divColumn1.add(createContentParagraphTwice(
                "The internal tumor barcode is " , failReport.sampleReport().tumorBarcode() , " and the internal blood barcode is "
                        , failReport.sampleReport().refBarcode()));
        divColumn1.add(createContentParagraph(
                "The tumor percentage estimated by Pathology UMC Utrecht is " , failReport.sampleReport().pathologyTumorPercentage()));
        return divColumn1;
    }

    @NotNull
    private Div createCPCTDRUPContentBodyColumn2() {
        Div divColumn2 = new Div();

        divColumn2.add(shallowSeqText());
        divColumn2.add(sampleArrivalDateText());
        divColumn2.add(recipientText());
        divColumn2.add(accreditationText());
        divColumn2.add(questionsText());
        return divColumn2;
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
        return createContentParagraph(
                "The tumor percentage based on molecular estimated is " , failReport.sampleReport().purityShallowSeq());
    }

    @NotNull
    private Paragraph sampleArrivalDateText() {
        return createContentParagraphTwice("The biopsies evaluated for this sample have arrived on ",
                DataUtil.formatDate(failReport.sampleReport().tumorArrivalDate()),
                " at ",
                ReportResources.HARTWIG_ADDRESS);
    }

    @NotNull
    private static Paragraph notSequencedText() {
        return createContentParagraph(
                "The received tumor biopsies were inadequate for whole genome sequencing. "
                        + "If available new tumor material can be provided for a new assessment.");
    }

    @NotNull
    private Paragraph recipientText() {
        return createContentParagraph(
                "This report is generated and verified by: " + failReport.user() + " and is addressed to " , failReport.sampleReport()
                        .addressee());
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
