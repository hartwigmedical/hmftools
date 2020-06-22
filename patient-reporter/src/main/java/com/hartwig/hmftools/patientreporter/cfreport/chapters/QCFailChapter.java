package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailType;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

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
        reportDocument.add(new Paragraph("The information regarding 'primary tumor location' and 'cancer subtype' is based on "
                + "information received from the originating hospital.").addStyle(ReportResources.subTextSmallStyle()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        LimsStudy study = LimsStudy.fromSampleId(failReport.sampleReport().tumorSampleId());

        switch (study) {
            case WIDE:
                reportDocument.add(createWIDEContentBody());
                break;
            case CORE:
                reportDocument.add(createCOREContentBody());
                break;
            case DRUP:
            case CPCT:
                reportDocument.add(createCPCTDRUPContentBody());
                break;
        }

        reportDocument.add(ReportSignature.createSignatureDiv(failReport.logoRVAPath(), failReport.signaturePath()).setMarginTop(15));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private static Div createFailReasonDiv(@NotNull QCFailReason failReason) {
        String reason = DataUtil.NA_STRING;
        String explanation = DataUtil.NA_STRING;
        String explanationDetail = DataUtil.NA_STRING;

        switch (failReason.type()) {
            case LOW_QUALITY_BIOPSY: {
                reason = "Insufficient quality of received biomaterial(s)";
                explanation = "The received biomaterial(s) did not meet the requirements that are needed for \n"
                        + "high quality Whole Genome Sequencing";
                break;
            }
            case TECHNICAL_FAILURE: {
                reason = "Technical failure";
                explanation = "Whole Genome Sequencing could not be successfully performed on the received biomaterial(s) \n"
                        + "due to technical problems";
            }
        }

        switch (failReason) {
            case INSUFFICIENT_DNA: {
                explanationDetail =
                        "The tumor percentage based on molecular estimation could not be determined due to insufficient tumor DNA";
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

    @NotNull
    private Table createWIDEContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createWIDESampleDetailsColumn()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createWIDEDisclaimerColumn()));
        return table;
    }

    @NotNull
    private Div createWIDESampleDetailsColumn() {
        Div div = createSampleDetailsDiv();
        div.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        div.add(reportIsBasedOnTumorSampleArrivedAt());
        div.add(reportIsBasedOnBloodSampleArrivedAt());
        div.add(resultsAreObtainedBetweenDates());
        if (failReport.sampleReport().hospitalPathologySampleId() != null) {
            div.add(reportIsForPathologySampleID());
        }
        if (failReport.reason().type() == QCFailType.LOW_QUALITY_BIOPSY) {
            div.add(sampleHasMolecularTumorPercentage());
        }

        return div;
    }

    @NotNull
    private Div createWIDEDisclaimerColumn() {
        Div div = createDisclaimerDiv();
        div.add(reportIsBasedOnBloodAndTumorSamples());
        div.add(testsArePerformedByAccreditedLab());
        div.add(reportIsVerifiedByAndAddressedTo());
        div.add(reportIsGeneratedByPatientReporterVersion());
        failReport.comments().ifPresent(comments -> div.add(createContentParagraphRed("Comments: " + comments)));
        div.add(resubmitSample());
        div.add(forQuestionsPleaseContactHMF());
        return div;
    }

    @NotNull
    private Table createCOREContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCORESampleDetailsColumn()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createCOREDisclaimerColumn()));
        return table;
    }

    @NotNull
    private Div createCORESampleDetailsColumn() {
        Div divColumn = createSampleDetailsDiv();
        divColumn.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        divColumn.add(reportIsBasedOnTumorSampleArrivedAt());
        divColumn.add(reportIsBasedOnBloodSampleArrivedAt());
        divColumn.add(resultsAreObtainedBetweenDates());
        divColumn.add(reportIsForHospitalPatientIDAndPathologySampleId());
        divColumn.add(reportIsForProjectAndSubmission());
        if (failReport.reason().type() == QCFailType.LOW_QUALITY_BIOPSY) {
            divColumn.add(sampleHasMolecularTumorPercentage());
        }
        return divColumn;
    }

    @NotNull
    private Div createCOREDisclaimerColumn() {
        Div divColumn = createDisclaimerDiv();
        divColumn.add(reportIsBasedOnBloodAndTumorSamples());
        divColumn.add(testsArePerformedByAccreditedLab());
        divColumn.add(reportIsVerifiedByAndAddressedTo());
        divColumn.add(reportIsGeneratedByPatientReporterVersion());
        failReport.comments().ifPresent(comments -> divColumn.add(createContentParagraphRed("Comments: " + comments)));
        divColumn.add(resubmitSample());
        divColumn.add(forQuestionsPleaseContactHMF());
        return divColumn;
    }

    @NotNull
    private Table createCPCTDRUPContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createCPCTDRUPSampleDetailsColumn()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createCPCTDRUPDisclaimerColumn()));
        return table;
    }

    @NotNull
    private Div createCPCTDRUPSampleDetailsColumn() {
        Div divColumn = createSampleDetailsDiv();
        divColumn.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        divColumn.add(reportIsBasedOnTumorSampleArrivedAt());
        divColumn.add(reportIsBasedOnBloodSampleArrivedAt());
        divColumn.add(resultsAreObtainedBetweenDates());
        if (failReport.reason().type() == QCFailType.LOW_QUALITY_BIOPSY) {
            divColumn.add(sampleHasMolecularTumorPercentage());
        }

        return divColumn;
    }

    @NotNull
    private Div createCPCTDRUPDisclaimerColumn() {
        Div divColumn = createDisclaimerDiv();
        divColumn.add(reportIsBasedOnBloodAndTumorSamples());
        divColumn.add(testsArePerformedByAccreditedLab());
        divColumn.add(reportIsVerifiedByAndAddressedTo());
        divColumn.add(reportIsGeneratedByPatientReporterVersion());
        failReport.comments().ifPresent(comments -> divColumn.add(createContentParagraphRed("Comments: " + comments)));
        divColumn.add(resubmitSample());
        divColumn.add(forQuestionsPleaseContactHMF());
        return divColumn;
    }

    @NotNull
    private Paragraph resubmitSample() {
        return createContentParagraph("If available new tumor/blood material can be provided for a new assessment, please contact ",
                "info@hartwigmedicalfoundation.nl");
    }

    @NotNull
    private Paragraph reportIsForHospitalPatientIDAndPathologySampleId() {
        if (failReport.sampleReport().hospitalPathologySampleId() != null) {
            return createContentParagraphTwice("The hospital patient ID is ",
                    failReport.sampleReport().hospitalPatientId(),
                    " and the pathology tissue ID is: ",
                    failReport.sampleReport().hospitalPathologySampleId());
        } else {
            return createContentParagraph("The hospital patient ID is ", failReport.sampleReport().hospitalPatientId());
        }
    }

    @NotNull
    private Paragraph reportIsForPathologySampleID() {
        return createContentParagraph("The pathology tissue ID is: ", failReport.sampleReport().hospitalPathologySampleId());
    }

    @NotNull
    private Paragraph reportIsForProjectAndSubmission() {
        return createContentParagraphTwice("The project name of the sample is ",
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
    private Paragraph reportIsBasedOnBloodAndTumorSamples() {
        return createContentParagraph("The results stated in this report are based on the tested tumor and blood sample.");
    }

    @NotNull
    private Paragraph reportIsGeneratedByPatientReporterVersion() {
        return createContentParagraphTwiceWithOneBold("This report is generated by patient reporter ",
                ReportResources.VERSION_REPORT,
                " based on ",
                failReport.qsFormNumber());
    }

    @NotNull
    private Paragraph forQuestionsPleaseContactHMF() {
        return createContentParagraph("For questions, please contact ", "info@hartwigmedicalfoundation.nl");
    }

    @NotNull
    private Paragraph sampleHasMolecularTumorPercentage() {
        String effectivePurity =
                failReport.wgsPurityString() != null ? failReport.wgsPurityString() : failReport.sampleReport().shallowSeqPurityString();
        if (effectivePurity.equals(Lims.PURITY_NOT_RELIABLE_STRING) || effectivePurity.equals(Lims.NOT_PERFORMED_STRING)) {
            return createContentParagraph("The tumor percentage based on molecular estimation", " could not be determined");
        } else {
            return createContentParagraph("The tumor percentage based on molecular estimation is ", effectivePurity);
        }
    }

    @NotNull
    private Paragraph samplesAreEvaluatedAtHMFAndWithSampleID() {
        return createContentParagraphTwice("The biomaterials are evaluated at ",
                ReportResources.HARTWIG_ADDRESS,
                " and are known under HMF sample ID  ",
                failReport.sampleReport().tumorSampleId());
    }

    @NotNull
    private Paragraph reportIsVerifiedByAndAddressedTo() {
        return createContentParagraph("This report is generated and verified by: " + failReport.user() + " and is addressed to ",
                failReport.sampleReport().addressee());
    }

    @NotNull
    private Paragraph testsArePerformedByAccreditedLab() {
        return createContentParagraph(
                "The results on this report are based on tests that are performed under ISO/ICE-17025:2005 accreditation.");
    }

    @NotNull
    private static Div createSampleDetailsDiv() {
        Div div = new Div();
        div.add(new Paragraph("Sample details").addStyle(ReportResources.smallBodyHeadingStyle()));
        return div;
    }

    @NotNull
    private static Div createDisclaimerDiv() {
        Div div = new Div();
        div.add(new Paragraph("Disclaimer").addStyle(ReportResources.smallBodyHeadingStyle()));
        return div;
    }


    @NotNull
    private static Paragraph createContentParagraphRed(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyleRed()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
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

    @NotNull
    private static Paragraph createContentParagraphTwiceWithOneBold(@NotNull String regularPart, @NotNull String boldPart,
            @NotNull String regularPart2, @NotNull String boldPart2) {
        return createContentParagraph(regularPart).add(new Text(boldPart).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(regularPart2)
                .add(new Text(boldPart2).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }
}
