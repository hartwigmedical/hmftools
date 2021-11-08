package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.itextpdf.io.IOException;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DetailsAndDisclaimerChapter implements ReportChapter {

    @NotNull
    private final AnalysedPatientReport patientReport;

    public DetailsAndDisclaimerChapter(@NotNull AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    public String name() {
        return "Sample details & disclaimers";
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public void render(@NotNull Document reportDocument) throws IOException {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createSampleDetailsDiv(patientReport)));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createDisclaimerDiv(patientReport)));
        reportDocument.add(table);

        reportDocument.add(ReportSignature.createSignatureDiv(patientReport.logoRVAPath(), patientReport.signaturePath()));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private static Div createSampleDetailsDiv(@NotNull AnalysedPatientReport patientReport) {
        SampleReport sampleReport = patientReport.sampleReport();
        LimsCohortConfig cohort = sampleReport.cohort();

        Div div = new Div();

        div.add(new Paragraph("Sample details").addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("The samples have been sequenced at ", ReportResources.HARTWIG_ADDRESS));
        div.add(createContentParagraph("The samples have been analyzed by Next Generation Sequencing using Whole Genome Sequencing"));

        div.add(generateHMFAndPathologySampleIDParagraph(patientReport.sampleReport()));

        div.add(generateGermlineChoicePatientParagraph(patientReport.sampleReport()));

        String earliestArrivalDate = sampleReport.earliestArrivalDate();
        div.add(createContentParagraphTwice("The results in this report have been obtained between ",
                DataUtil.formatNullableString(earliestArrivalDate),
                " and ",
                ReportResources.REPORT_DATE));

        div.add(createContentParagraphTwice("This experiment is performed on the tumor sample which arrived on ",
                DataUtil.formatDate(sampleReport.tumorArrivalDate()),
                " with internal tumor barcode ",
                sampleReport.tumorSampleBarcode()));
        div.add(createContentParagraphTwice("This experiment is performed on the blood sample which arrived on ",
                DataUtil.formatDate(sampleReport.refArrivalDate()),
                " with internal blood barcode ",
                DataUtil.formatNullableString(sampleReport.refSampleBarcode())));
        div.add(createContentParagraph("The results stated in this report are based on the tested tumor and blood sample."));
        div.add(createContentParagraph("This experiment is performed according to lab procedures: ", sampleReport.labProcedures()));
        String whoVerified = "This report was generated " + patientReport.user();

        div.add(createContentParagraph(whoVerified));
        div.add(createContentParagraph("This report is addressed to: ", sampleReport.addressee()));

        if (cohort.requireHospitalId()) {
            div.add(createContentParagraph("The hospital patient ID is: ", sampleReport.hospitalPatientId()));
        }

        if (cohort.requireHospitalPersonsRequester()) {
            div.add(createContentParagraphTwice("The project name of sample is: ",
                    sampleReport.projectName(),
                    " and the submission ID is ",
                    sampleReport.submissionId()));
        }

        patientReport.comments().ifPresent(comments -> div.add(createContentParagraphRed("Comments: " + comments)));

        return div;
    }

    @NotNull
    private static Div createDisclaimerDiv(@NotNull AnalysedPatientReport patientReport) {
        String pipelineVersion = patientReport.pipelineVersion() == null ? "No pipeline version is known" : patientReport.pipelineVersion();
        Div div = new Div();

        div.add(new Paragraph("Disclaimer").addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("The data on which this report is based is generated "
                + "from tests that are performed under ISO/ICE-17025:2017 TESTING L633 accreditation and have passed all internal quality controls."));
        div.add(createContentParagraphTwice("This report is generated by patient reporter ",
                ReportResources.VERSION_REPORT,
                " based on ",
                patientReport.qsFormNumber() + "."));
        div.add(createContentParagraph("UDI-DI: ", patientReport.udiDi() + "."));
        div.add(createContentDivWithLink("The OncoAct user manual can be found at ", ReportResources.MANUAL + ".", ReportResources.MANUAL + "."));
        div.add(createContentParagraph("This report is based on pipeline version ", pipelineVersion + "."));
        div.add(createContentParagraph("The ‘primary tumor location’ and ‘primary tumor type’ have influence on the "
                + "clinical evidence/study matching. No check is performed to verify the received information."));
        div.add(createContentParagraph("The conclusion of this report is based solely on the results of the DNA sequencing of the tumor "
                + "and the received tumor type. Final interpretation of the clinical consequence of this report should therefore "
                + "always be performed by the treating physician."));
        div.add(createContentParagraph("Based on a tumor purity of at least 20%, the test has a sensitivity of >95% for detection of "
                + "somatic variants and >95% for detection of translocations and gene copy number changes."));
        div.add(createContentParagraph("For feedback or complaints please contact ", ReportResources.CONTACT_EMAIL_QA + "."));
        div.add(createContentParagraph("For questions about the contents of this report, please contact ",
                ReportResources.CONTACT_EMAIL_GENERAL + "."));

        return div;
    }

    @NotNull
    private static Div createContentDivWithLink(@NotNull String string1, @NotNull String string2, @NotNull String link) {
        Div div = new Div();

        div.add(createParaGraphWithLink(string1, string2, link));
        return div;
    }

    @NotNull
    private static Paragraph createParaGraphWithLink(@NotNull String string1, @NotNull String string2, @NotNull String link) {
        return new Paragraph(string1).addStyle(ReportResources.subTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string2).addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(link)))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph generateHMFAndPathologySampleIDParagraph(@NotNull SampleReport sampleReport) {
        if (sampleReport.hospitalPathologySampleId() != null && sampleReport.cohort().requireHospitalPAId()) {
            return createContentParagraphTwice("The HMF sample ID is: ",
                    sampleReport.tumorSampleId(),
                    " and the pathology tissue ID is: ",
                    sampleReport.hospitalPathologySampleId());
        } else {
            return createContentParagraph("The HMF sample ID is: ", sampleReport.tumorSampleId());
        }
    }

    @NotNull
    private static Paragraph generateGermlineChoicePatientParagraph(@NotNull SampleReport sampleReport) {
        return createContentParagraph("The germline reporting choice of this patient is: ",
                sampleReport.germlineReportingLevel().display());
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

}
