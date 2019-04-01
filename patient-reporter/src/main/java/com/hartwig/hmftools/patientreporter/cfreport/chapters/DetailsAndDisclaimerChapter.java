package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.DataUtility;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.hartwig.hmftools.patientreporter.report.pages.SampleDetailsPage;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.*;
import com.itextpdf.layout.property.UnitValue;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.net.MalformedURLException;

public class DetailsAndDisclaimerChapter extends ReportChapter {

    private static final Logger LOGGER = LogManager.getLogger(SampleDetailsPage.class);

    @Override
    public String getName() {
        return "Sample details & disclaimers";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ClosingChapter;
    }

    @Override
    protected void renderChapterContent(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) throws IOException {

        // Main content
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 0.1f, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(createSampleDetailsDiv(patientReport)));
        table.addCell(TableHelper.getLayoutCell()); // Spacer
        table.addCell(TableHelper.getLayoutCell()
                .add(createDisclaimerDiv()));
        reportDocument.add(table);

        // End of report text
        reportDocument.add(new Paragraph("— End of report —")
                .setMarginTop(50)
                .addStyle(ReportResources.smallBodyTextStyle()));

        reportDocument.add(createSignatureDiv(patientReport.logoRVAPath(), patientReport.signaturePath()).setPaddingTop(80));

    }

    @NotNull
    private static Div createSampleDetailsDiv(@NotNull final AnalysedPatientReport patientReport) {

        final SampleReport sampleReport = patientReport.sampleReport();

        String recipient = sampleReport.recipient();
        if (recipient == null) {
            LOGGER.warn("No recipient address present for sample " + sampleReport.sampleId());
            recipient = DataUtility.NAString;
        }

        Div div = new Div();

        // Heading
        div.add(new Paragraph("Sample details")
                .addStyle(ReportResources.smallBodyHeadingStyle()));

        // Content
        div.add(createContentParagraph("The samples have been sequenced at ", ReportResources.HARTWIG_ADDRESS));
        div.add(createContentParagraph("The samples have been analyzed by Next Generation Sequencing "));
        div.add(createContentParagraph("This experiment is performed on the tumor sample which arrived on ", DataUtility.formatDate(sampleReport.tumorArrivalDate())));
        div.add(createContentParagraph("The pathology tumor percentage for this sample is " + sampleReport.pathologyTumorPercentage()));
        div.add(createContentParagraph("This experiment is performed on the blood sample which arrived on ", DataUtility.formatDate(sampleReport.bloodArrivalDate())));
        div.add(createContentParagraph("This experiment is performed according to lab procedures: " + sampleReport.labProcedures()));
        div.add(createContentParagraph("This report is generated and verified by: " + patientReport.user()));
        div.add(createContentParagraph("This report is addressed at: " + recipient));
        patientReport.comments().ifPresent(comments -> div.add(createContentParagraph("Comments: " + comments)));

        return div;

    }

    @NotNull
    private static Div createDisclaimerDiv() {

        Div div = new Div();

        // Heading
        div.add(new Paragraph("Disclaimer")
                .addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("The data on which this report is based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation."));
        div.add(createContentParagraph("The analysis done for this report has passed all internal quality controls."));
        div.add(createContentParagraph("For feedback or complaints please contact ", ReportResources.CONTACT_EMAIL_QA));
        div.add(createContentParagraph("For general questions, please contact us at ", ReportResources.CONTACT_EMAIL_GENERAL));

        return div;

    }

    @NotNull
    private static Div createSignatureDiv(@NotNull String rvaLogoPath, @NotNull String signaturePath) throws IOException {

        Div div = new Div();
        div.setKeepTogether(true);

        // Add RVA logo
        try {
            final Image rvaLogo = new Image(ImageDataFactory.create(rvaLogoPath));
            rvaLogo.setMaxHeight(58);
            if (rvaLogo != null) {
                div.add(rvaLogo);
            }
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read RVA logo image at " + rvaLogoPath);
        }

        // Add signature text
        Paragraph signatureText = new Paragraph()
                .setFont(ReportResources.getFontBold())
                .setFontSize(10)
                .setFontColor(ReportResources.PALETTE_BLACK);

        signatureText.add(ReportResources.SIGNATURE_NAME + ",\n");
        signatureText.add(new Text(ReportResources.SIGNATURE_TITLE).setFont(ReportResources.getFontRegular()));
        div.add(signatureText);

        // Add signature image
        try {
            final Image signatureImage = new Image(ImageDataFactory.create(signaturePath));
            signatureImage.setMaxHeight(60);
            signatureImage.setMarginTop(-15); // Set negative margin so the signature slightly overlaps the signature text
            signatureImage.setMarginLeft(10);
            div.add(signatureImage);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read signature image at " + signaturePath);
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
    private static Paragraph createContentParagraph(@NotNull String regularPart, @NotNull String boldPart) {
        return createContentParagraph(regularPart)
                .add(new Text(boldPart)
                        .addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

}
