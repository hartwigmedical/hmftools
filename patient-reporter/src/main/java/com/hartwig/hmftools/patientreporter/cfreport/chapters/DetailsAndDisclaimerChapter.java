package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.itextpdf.io.IOException;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class DetailsAndDisclaimerChapter implements ReportChapter {

    @NotNull
    private final AnalysedPatientReport patientReport;

    public DetailsAndDisclaimerChapter(@NotNull AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
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
        table.addCell(TableUtil.createLayoutCell().add(createDisclaimerDiv()));
        reportDocument.add(table);

        reportDocument.add(ReportSignature.createSignatureDiv(patientReport.logoRVAPath(), patientReport.signaturePath()));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private static Div createSampleDetailsDiv(@NotNull AnalysedPatientReport patientReport) {
        SampleReport sampleReport = patientReport.sampleReport();
        LimsStudy study = LimsStudy.fromSampleId(patientReport.sampleReport().tumorSampleId());

        Paragraph sampleIdentificationLineOnReport;
        if (patientReport.sampleReport().hospitalPathologySampleId() != null) {
            sampleIdentificationLineOnReport = createContentParagraphTwice("The HMF sample ID is: ",
                    patientReport.sampleReport().tumorSampleId(),
                    " and the pathology tissue ID is: ",
                    patientReport.sampleReport().hospitalPathologySampleId());
        } else {
            sampleIdentificationLineOnReport =
                    createContentParagraph("The HMF sample ID is: ", patientReport.sampleReport().tumorSampleId());
        }

        Div div = new Div();

        div.add(new Paragraph("Sample details").addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("The samples have been sequenced at ", ReportResources.HARTWIG_ADDRESS));
        div.add(createContentParagraph("The samples have been analyzed by Next Generation Sequencing using Whole Genome Sequencing"));

        String earliestArrivalDate = sampleReport.earliestArrivalDate();
        div.add(createContentParagraphTwice("The results in this report have been obtained between ",
                earliestArrivalDate != null ? earliestArrivalDate : DataUtil.NA_STRING,
                " and ",
                ReportResources.REPORT_DATE));

        div.add(sampleIdentificationLineOnReport);
        div.add(createContentParagraphTwice("This experiment is performed on the tumor sample which arrived on ",
                DataUtil.formatDate(sampleReport.tumorArrivalDate()),
                " with internal tumor barcode ",
                sampleReport.tumorSampleBarcode()));
        div.add(createContentParagraphTwice("This experiment is performed on the blood sample which arrived on ",
                DataUtil.formatDate(sampleReport.refArrivalDate()),
                " with internal blood barcode ",
                sampleReport.refSampleBarcode()));
        div.add(createContentParagraph("This experiment is performed according to lab procedures: ", sampleReport.labProcedures()));
        String whoVerified = "This report is generated and verified by: " + patientReport.user();
        if (!patientReport.clinicalSummary().isEmpty()) {
            whoVerified += " and Paul Roepman";
        }
        div.add(createContentParagraph(whoVerified));
        div.add(createContentParagraph("This report is addressed to: ", sampleReport.addressee()));

        if (study == LimsStudy.CORE) {
            div.add(createContentParagraph("The hospital patient ID is: ", sampleReport.hospitalPatientId()));
            div.add(createContentParagraphTwice("The project name of sample is: ",
                    sampleReport.projectName(),
                    " and the submission ID is ",
                    sampleReport.submissionId()));
            div.add(createContentParagraphRequest(sampleReport));
        }
        patientReport.comments().ifPresent(comments -> div.add(createContentParagraph("Comments: " + comments)));

        return div;
    }

    @NotNull
    private static Div createDisclaimerDiv() {
        Div div = new Div();

        div.add(new Paragraph("Disclaimer").addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("This report is generated by patient reporter ", ReportResources.VERSION_REPORT));
        div.add(createContentParagraph("The data on which this report is based is generated "
                + "from tests that are performed under ISO/ICE-17025:2005 accreditation and have passed all internal quality controls."));
        div.add(createContentParagraph("The results stated in this report are based on the tested tumor and blood sample."));
        div.add(createContentParagraph("The ‘primary tumor location’ and ‘cancer subtype’ have influence on the "
                + "clinical evidence/study matching. No check is performed to verify the received information."));
        div.add(createContentParagraph("The conclusion of this report is based solely on the results of the DNA sequencing of the tumor"
                + " and the received tumor type. All other patient/tumor characteristics that might influence the interpretation of "
                + "these results, are not considered. Final interpretation of the clinical consequence of this report should therefore "
                + "always be performed by the treating physician."));
        div.add(createContentParagraph("Based on a tumor purity of at least 30%, the test has a sensitivity of >95% for detection "
                + "of somatic variants and >95% for detection of translocations and gene copy number changes. For samples "
                + "with a purity above 20%, the test has a sensitivity of >90%."));
        div.add(createContentParagraphTwice("For feedback or complaints please contact ",
                ReportResources.CONTACT_EMAIL_QA,
                " and for general questions, please contact ",
                ReportResources.CONTACT_EMAIL_GENERAL));

        return div;
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
    private static Paragraph createContentParagraphRequest(@NotNull SampleReport sampleReport) {
        String requesterName = sampleReport.hospitalContactData().requesterName();
        String requesterEmail = sampleReport.hospitalContactData().requesterEmail();

        return createContentParagraph("The requester is: ").add(new Text(requesterName).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(new Text(" (" + requesterEmail + ")").addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }
}
