package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;
import org.jetbrains.annotations.NotNull;

public class DetailsAndDisclaimerChapter extends ReportChapter {

    private final static Style BODY_TEXT_HEADER_STYLE = ReportResources.smallBodyHeadingStyle();
    private final static Style BODY_TEXT_STYLE = ReportResources.smallBodyTextStyle();
    private final static Style BODY_TEXT_BOLD_STYLE = ReportResources.smallBodyBoldTextStyle();


    @Override
    public String getName() {
        return "Sample details & disclaimers";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ClosingChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {

        // Main content
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 0.1f, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(getSampleDetailsDiv()));
        table.addCell(TableHelper.getLayoutCell()); // Spacer
        table.addCell(TableHelper.getLayoutCell()
                .add(getDisclaimerDiv()));
        report.add(table);

        // End of report text
        report.add(new Paragraph("— End of report —")
                .setMarginTop(50)
                .addStyle(BODY_TEXT_STYLE));

        // @TODO Add signature


    }

    private static final Div getSampleDetailsDiv() {

        Div div = new Div();

        // Heading
        div.add(new Paragraph("Sample details")
                .addStyle(BODY_TEXT_HEADER_STYLE));

        // Content
        div.add(getContentParagraph("The samples have been sequenced at ", ReportResources.HARTWIG_ADDRESS));
        div.add(getContentParagraph("The samples have been analyzed by Next Generation Sequencing "));
        div.add(getContentParagraph("This experiment is performed on the tumor sample which arrived on ", "05-Jan-2018"));
        div.add(getContentParagraph("The pathology tumor percentage for this sample is 80%"));
        div.add(getContentParagraph("This experiment is performed on the blood sample which arrived on ", "01-Jan-2018"));
        div.add(getContentParagraph("This experiment is performed according to lab procedures: PREP013V23-QC037V20-SEQ008V25"));
        div.add(getContentParagraph("This report is generated and verified by: korneelduyvesteyn"));
        div.add(getContentParagraph("This report is addressed at: HMF Testing Center"));
        div.add(getContentParagraph("Comments: this is a test report and is based off COLO829"));

        return div;

    }

    private static final Div getDisclaimerDiv() {

        Div div = new Div();

        // Heading
        div.add(new Paragraph("Disclaimer")
                .addStyle(BODY_TEXT_HEADER_STYLE));

        div.add(getContentParagraph("The data on which this report is based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation."));
        div.add(getContentParagraph("The analysis done for this report has passed all internal quality controls."));
        div.add(getContentParagraph("For feedback or complaints please contact ", ReportResources.CONTACT_EMAIL_QA));
        div.add(getContentParagraph("For general questions, please contact us at ", ReportResources.CONTACT_EMAIL_GENERAL));

        return div;

    }

    @NotNull
    private final static Paragraph getContentParagraph(@NotNull String text) {
        return new Paragraph(text)
                .addStyle(BODY_TEXT_STYLE)
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private final static Paragraph getContentParagraph(@NotNull String regularPart, @NotNull String boldPart) {

        return getContentParagraph(regularPart)
                .add(new Text(boldPart)
                        .addStyle(BODY_TEXT_BOLD_STYLE))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

}
