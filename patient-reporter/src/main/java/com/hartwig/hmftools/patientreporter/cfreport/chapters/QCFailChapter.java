package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsCoreCohort;
import com.hartwig.hmftools.common.lims.LimsStudy;
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
        return failReport.isCorrectedReport() ? failReport.reason().title() + " (Corrected)" : failReport.reason().title();
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
        if (failReport.sampleReport().hospitalQuery().hospitalPIWithAddress() == null) {
            throw new IllegalStateException("No recipient address present for sample " + failReport.sampleReport().tumorSampleId());
        }

        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(failReport.sampleReport().primaryTumorLocationString(),
                failReport.sampleReport().cancerSubTypeString(),
                contentWidth()));
        reportDocument.add(new Paragraph("The information regarding 'primary tumor location' and 'cancer subtype' is based on "
                + "information received from the originating hospital.").addStyle(ReportResources.subTextSmallStyle()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        LimsStudy study = LimsStudy.fromSampleId(failReport.sampleReport().tumorSampleId());

        switch (study) {
            case CORE:
                reportDocument.add(createCOREContentBody());
                break;
            case WIDE:
                reportDocument.add(createWIDEContentBody());
                break;
            default:
                reportDocument.add(createCPCTDRUPContentBody());
        }

        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setMarginTop(15));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason) {
        String title;
        String reason;
        String explanation;

        switch (failReason) {
            case LOW_DNA_YIELD: {
                title = "Notification failed sample";
                reason = "Insufficient biopsy/tissue quality";
                explanation = "The received biopsy/tissue sample did not meet the requirements that are needed for \n high quality "
                        + "Whole Genome Sequencing";
                break;
            }
            case POST_ANALYSIS_FAIL: {
                title = "Notification failed sample";
                reason = "Insufficient biopsy/tissue quality";
                explanation = "The received biopsy/tissue sample did not meet the requirements that are needed for \n high quality "
                        + "Whole Genome Sequencing";
                break;
            }
            case SHALLOW_SEQ_LOW_PURITY: {
                title = "Notification failed sample";
                reason = "Insufficient biopsy/tissue quality";
                explanation = "The received biopsy/tissue sample did not meet the requirements that are needed for \n high quality "
                        + "Whole Genome Sequencing";
                break;
            }
            case INSUFFICIENT_TISSUE: {
                title = "Notification failed sample";
                reason = "Insufficient biopsy/tissue quality";
                explanation = "The received biopsy/tissue sample did not meet the requirements that are needed for \n high quality "
                        + "Whole Genome Sequencing";
                break;
            }
            case BELOW_DETECTION_THRESHOLD: {
                title = "Notification failed sample";
                reason = "Insufficient biopsy/tissue quality";
                explanation = "The received biopsy/tissue sample did not meet the requirements that are needed for \n high quality "
                        + "Whole Genome Sequencing";
                break;
            }
            case LAB_FAILURE: {
                title = "Notification failed sample";
                reason = "Technical failure";
                explanation = "Whole Genome Sequencing could not be successfully performed on the received biopsy \n "
                        + "due to technical problems";
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
        table.addCell(TableUtil.createLayoutCell().add(createWIDEContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createWIDEContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createWIDEContentBodyColumn1() {
        Div divColumn = new Div();
        divColumn.add(sampleNotAdequateForWGS());
        divColumn.add(resubmitInSameStudyWithSameNumber());
        if (failReport.sampleReport().hospitalPathologySampleId() != null) {
            divColumn.add(reportIsForPathologyTissueID());
        }
        divColumn.add(resultsAreObtainedBetweenDates());
        divColumn.add(reportIsBasedOnTumorSampleArrivedAt());
        divColumn.add(reportIsBasedOnBloodSampleArrivedAt());
        divColumn.add(sampleHasMolecularTumorPercentage());
        failReport.comments().ifPresent(comments -> divColumn.add(createContentParagraph("Comments: " + comments)));

        return divColumn;
    }

    @NotNull
    private Div createWIDEContentBodyColumn2() {
        Div divColumn = new Div();
        divColumn.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        divColumn.add(reportIsVerifiedByAndAddressedAt());
        divColumn.add(reportIsBasedOnBloodAndTumorSamples());
        divColumn.add(reportIsGeneratedByPatientReporterVersion());
        divColumn.add(testsArePerformedByAccreditedLab());
        divColumn.add(forQuestionsPleaseContactHMF());
        return divColumn;
    }

    @NotNull
    private Table createCOREContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCOREContentBodyColumn1()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createCOREContentBodyColumn2()));
        return table;
    }

    @NotNull
    private Div createCOREContentBodyColumn1() {
        LimsCoreCohort coreCohort = LimsCoreCohort.fromSampleId(failReport.sampleReport().tumorSampleId());
        Div divColumn = new Div();
        divColumn.add(sampleNotAdequateForWGS());
        divColumn.add(resubmitInSameDVOWithProjectName());
        if (failReport.sampleReport().hospitalPathologySampleId() != null) {
            divColumn.add(reportIsForPathologyTissueID());
        }
        divColumn.add(reportIsForHospitalPatientID());
        divColumn.add(reportIsForProjectAndSubmission());
        divColumn.add(resultsAreObtainedBetweenDates());
        divColumn.add(reportIsBasedOnTumorSampleArrivedAt());
        divColumn.add(reportIsBasedOnBloodSampleArrivedAt());
        divColumn.add(sampleHasMolecularTumorPercentage());
        failReport.comments().ifPresent(comments -> divColumn.add(createContentParagraph("Comments: " + comments)));
        return divColumn;
    }

    @NotNull
    private Div createCOREContentBodyColumn2() {
        Div divColumn = new Div();
        divColumn.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        divColumn.add(reportIsVerifiedByAndAddressedAt());
        divColumn.add(reportIsBasedOnBloodAndTumorSamples());
        divColumn.add(reportIsRequestedBy());
        divColumn.add(reportIsGeneratedByPatientReporterVersion());
        divColumn.add(testsArePerformedByAccreditedLab());
        divColumn.add(forQuestionsPleaseContactHMF());
        return divColumn;
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
        Div divColumn = new Div();
        divColumn.add(sampleNotAdequateForWGS());
        divColumn.add(resubmitInSameStudyWithSameNumber());
        divColumn.add(resultsAreObtainedBetweenDates());
        divColumn.add(reportIsBasedOnTumorSampleArrivedAt());
        divColumn.add(reportIsBasedOnBloodSampleArrivedAt());
        divColumn.add(sampleHasMolecularTumorPercentage());
        failReport.comments().ifPresent(comments -> divColumn.add(createContentParagraph("Comments: " + comments)));

        return divColumn;
    }

    @NotNull
    private Div createCPCTDRUPContentBodyColumn2() {
        Div divColumn = new Div();
        divColumn.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        divColumn.add(reportIsVerifiedByAndAddressedAt());
        divColumn.add(reportIsBasedOnBloodAndTumorSamples());
        divColumn.add(reportIsGeneratedByPatientReporterVersion());
        divColumn.add(testsArePerformedByAccreditedLab());
        divColumn.add(forQuestionsPleaseContactHMF());
        return divColumn;
    }

    @NotNull
    private Paragraph sampleNotAdequateForWGS() {
        return createContentParagraph("The received tumor biopsies were inadequate for whole genome sequencing. ");
    }

    @NotNull
    private Paragraph resubmitInSameStudyWithSameNumber() {
        return createContentParagraph(
                "If available new tumor material can be provided for a new assessment, please resubmit using the same " + failReport.study()
                        .studyName() + "-number. " + "If additional material cannot be provided the patient will not be "
                        + "evaluable for the " + failReport.study().studyCode() + " study.");
    }

    @NotNull
    private Paragraph resubmitInSameDVOWithProjectName() {
        return createContentParagraph("If available new tumor material can be provided for a new assessment, "
                + "please resubmit using the same DVO with project name " + failReport.sampleReport().projectName() + ".");
    }

    @NotNull
    private Paragraph reportIsForHospitalPatientID() {
        return createContentParagraph("The hospital patient ID is ", failReport.sampleReport().hospitalPatientId());
    }

    @NotNull
    private Paragraph reportIsForPathologyTissueID() {
        return createContentParagraph("The tissue ID of pathology is: ", failReport.sampleReport().hospitalPathologySampleId());
    }

    @NotNull
    private Paragraph reportIsForProjectAndSubmission() {
        return createContentParagraphTwice("The project name of sample is ",
                failReport.sampleReport().projectName(),
                " and the submission ID is ",
                failReport.sampleReport().submissionId());
    }

    @NotNull
    private Paragraph resultsAreObtainedBetweenDates() {
        String earliestArrivalDate = failReport.sampleReport().earliestArrivalDate();
        return createContentParagraphTwice("The results in this report have been obtained between ",
                earliestArrivalDate != null ? earliestArrivalDate : DataUtil.NA_STRING,
                " and ",
                ReportResources.REPORT_DATE);
    }

    @NotNull
    private Paragraph reportIsBasedOnTumorSampleArrivedAt() {
        return createContentParagraphTwice("This experiment is performed on the tumor sample which arrived on ",
                DataUtil.formatDate(failReport.sampleReport().tumorArrivalDate()),
                " with internal tumor barcode ",
                failReport.sampleReport().tumorSampleBarcode());
    }

    @NotNull
    private Paragraph reportIsBasedOnBloodSampleArrivedAt() {
        return createContentParagraphTwice("This experiment is performed on the blood sample which arrived on ",
                DataUtil.formatDate(failReport.sampleReport().refArrivalDate()),
                " with internal blood barcode ",
                failReport.sampleReport().refSampleBarcode());
    }

    @NotNull
    private Paragraph reportIsRequestedBy() {
        String requesterName = failReport.sampleReport().hospitalQuery().requestName();
        String requesterEmail = failReport.sampleReport().hospitalQuery().requestEmail();
        return createContentParagraph("The requester is : ").add(new Text(requesterName).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(new Text(" (" + requesterEmail + ")").addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private Paragraph reportIsBasedOnBloodAndTumorSamples() {
        return createContentParagraph("The results stated in this report are based on the tested tumor and blood sample.");
    }

    @NotNull
    private Paragraph reportIsGeneratedByPatientReporterVersion() {
        return createContentParagraph("This report is based by patient reporter ", ReportResources.VERSION_REPORT);
    }

    @NotNull
    private Paragraph forQuestionsPleaseContactHMF() {
        return createContentParagraph("For questions, please contact ", "info@hartwigmedicalfoundation.nl");
    }

    @NotNull
    private Paragraph sampleHasMolecularTumorPercentage() {
        return createContentParagraph("The tumor percentage based on molecular estimation is ",
                failReport.sampleReport().purityShallowSeq());
    }

    @NotNull
    private Paragraph samplesAreEvaluatedAtHMFAndWithSampleID() {
        return createContentParagraphTwice("The biopsies are evaluated at ",
                ReportResources.HARTWIG_ADDRESS,
                " and are known under HMF sample ID  ",
                failReport.sampleReport().tumorSampleId());
    }

    @NotNull
    private Paragraph reportIsVerifiedByAndAddressedAt() {
        String addressee = failReport.sampleReport().hospitalQuery().hospitalPIWithAddress();
        assert addressee != null; // Has been checked prior to calling this function.
        return createContentParagraph("This report is generated and verified by: " + failReport.user() + " and is addressed to ",
                addressee);
    }

    @NotNull
    private Paragraph testsArePerformedByAccreditedLab() {
        return createContentParagraph(
                "The results on this report are based on tests that are performed under ISO/ICE-17025:2005 accreditation.");
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
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
