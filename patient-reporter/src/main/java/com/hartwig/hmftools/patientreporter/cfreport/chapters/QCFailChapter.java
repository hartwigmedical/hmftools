package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.common.utils.DataUtil;
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
        reportDocument.add(TumorLocationAndTypeTable.createBiopsyLocationAndTumorLocation(failReport.sampleReport()
                .primaryTumorLocationString(), failReport.sampleReport().biopsyLocation(), contentWidth()));
        reportDocument.add(new Paragraph());
        reportDocument.add(TumorLocationAndTypeTable.createTumorType(failReport.sampleReport().primaryTumorTypeString(),
                contentWidth()));

        reportDocument.add(new Paragraph("The information regarding 'primary tumor location', 'primary tumor type' 'biopsy location'"
                + " is based on information received from the originating hospital.").addStyle(ReportResources.subTextSmallStyle()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createFailReasonDiv(failReport.reason()));
        reportDocument.add(LineDivider.createLineDivider(contentWidth()));

        reportDocument.add(createContentBody());
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
                        "The tumor percentage based on molecular estimation could not be determined due to insufficient tumor DNA.";
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
    private Table createContentBody() {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createSampleDetailsColumn()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createDisclaimerColumn()));
        return table;
    }

    @NotNull
    private Div createSampleDetailsColumn() {
        Div div = createSampleDetailsDiv();
        div.add(samplesAreEvaluatedAtHMFAndWithSampleID());
        div.add(reportIsBasedOnTumorSampleArrivedAt());
        div.add(reportIsBasedOnBloodSampleArrivedAt());
        div.add(resultsAreObtainedBetweenDates());
        if (failReport.sampleReport().cohort().requireHospitalPAId() && !failReport.sampleReport().cohort().requireHospitalId()) {
            if (failReport.sampleReport().hospitalPathologySampleId() != null){
                div.add(reportIsForPathologySampleID());
            }
        }
        if (failReport.sampleReport().cohort().requireHospitalPAId() && failReport.sampleReport().cohort().requireHospitalId()) {
            if (failReport.sampleReport().hospitalPathologySampleId() != null &&  failReport.sampleReport().hospitalPatientId() != null) {
                div.add(reportHospitalPatientIDAndPathologySampleId());
            }
        }

        if (failReport.sampleReport().cohort().requireHospitalPersonsRequester()) {
            div.add(reportIsForProjectAndSubmission());
        }

        if (failReport.reason().type() == QCFailType.LOW_QUALITY_BIOPSY) {
            div.add(sampleHasMolecularTumorPercentage());
        }

        return div;
    }

    @NotNull
    private Div createDisclaimerColumn() {
        Div div = createDisclaimerDiv();
        div.add(reportIsBasedOnBloodAndTumorSamples());
        div.add(testsArePerformedByAccreditedLab());
        div.add(testsArePerformedUnderUNI());
        div.add(testsManual());
        div.add(reportIsVerifiedByAndAddressedTo());
        div.add(reportIsGeneratedByPatientReporterVersion());
        failReport.comments().ifPresent(comments -> div.add(createContentParagraphRed("Comments: " + comments)));
        div.add(resubmitSample());
        div.add(forQuestionsPleaseContactHMF());
        return div;
    }

    @NotNull
    private Paragraph resubmitSample() {
        return createContentParagraph("If available new biomaterial(s) can be provided for a new assessment, please contact ",
                "info@hartwigmedicalfoundation.nl");
    }

    @NotNull
    private Paragraph reportHospitalPatientIDAndPathologySampleId() {
        return createContentParagraphTwice("The hospital patient ID is ",
                failReport.sampleReport().hospitalPatientId(),
                " and the pathology tissue ID is: ",
                failReport.sampleReport().hospitalPathologySampleId());

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
                DataUtil.formatNullableString(failReport.sampleReport().refSampleBarcode()));
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
                failReport.qsFormNumber() + ".");
    }

    @NotNull
    private Paragraph forQuestionsPleaseContactHMF() {
        return createContentParagraph("For questions regarding the results described in this report, please contact ",
                ReportResources.CONTACT_EMAIL_GENERAL);
    }

    @NotNull
    private Paragraph sampleHasMolecularTumorPercentage() {
        String effectivePurity =
                failReport.wgsPurityString() != null ? failReport.wgsPurityString() : failReport.sampleReport().shallowSeqPurityString();
        if (effectivePurity.equals(Lims.PURITY_NOT_RELIABLE_STRING) || effectivePurity.equals(Lims.NOT_PERFORMED_STRING)) {
            return createContentParagraph("The tumor percentage based on molecular estimation", " could not be determined.");
        } else {
            return createContentParagraph("The tumor percentage based on molecular estimation is ", effectivePurity);
        }
    }

    @NotNull
    private Paragraph samplesAreEvaluatedAtHMFAndWithSampleID() {
        return createContentParagraphTwice("The biomaterials are evaluated at ",
                ReportResources.HARTWIG_ADDRESS,
                " and are known as HMF sample ID  ",
                failReport.sampleReport().tumorSampleId());
    }

    @NotNull
    private Paragraph reportIsVerifiedByAndAddressedTo() {
        return createContentParagraph("This report is generated and verified by: " + failReport.user() + " and is addressed to ",
                failReport.sampleReport().addressee() + ".");
    }

    @NotNull
    private Paragraph testsArePerformedByAccreditedLab() {
        return createContentParagraph(
                "The results on this report are based on tests that are performed under ISO/ICE-17025:2017 TESTING L633 accreditation.");
    }

    @NotNull
    private Paragraph testsArePerformedUnderUNI() {
        return createContentParagraph("UDI-DI: ", ReportResources.UDI_DI + ".");
    }

    @NotNull
    private Paragraph testsManual() {
        return createContentParagraph("The OncoAct user manual can be found at ", ReportResources.MANUAL + ".");
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
