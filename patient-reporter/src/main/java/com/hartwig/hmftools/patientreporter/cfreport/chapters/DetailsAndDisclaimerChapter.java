package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.BodyText;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.*;
import com.itextpdf.layout.property.UnitValue;

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
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(getSampleDetailsDiv()));
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
        div.add(new Paragraph("The samples have been sequenced at ")
                .addStyle(BODY_TEXT_STYLE)
                .add(new Text(ReportResources.HARTWIG_ADDRESS)
                        .addStyle(BODY_TEXT_BOLD_STYLE)));

        div.add(new Paragraph("The samples have been analyzed by Next Generation Sequencing")
                        .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("This experiment is performed on the tumor sample which arrived on ")
                .addStyle(BODY_TEXT_STYLE)
                .add(new Text("05-Jan-2018")
                        .addStyle(BODY_TEXT_BOLD_STYLE)));

        div.add(new Paragraph("The pathology tumor percentage for this sample is 80%")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("This experiment is performed on the blood sample which arrived on ")
                .addStyle(BODY_TEXT_STYLE)
                .add(new Text("01-Jan-2018")
                        .addStyle(BODY_TEXT_BOLD_STYLE)));

        div.add(new Paragraph("This experiment is performed according to lab procedures: PREP013V23-QC037V20-SEQ008V25")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("This report is generated and verified by: korneelduyvesteyn")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("This report is addressed at: HMF Testing Center")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("Comments: this is a test report and is based off COLO829")
                .addStyle(BODY_TEXT_STYLE));

        return div;

    }

    private static final Div getDisclaimerDiv() {

        Div div = new Div();

        // Heading
        div.add(new Paragraph("Disclaimer")
                .addStyle(BODY_TEXT_HEADER_STYLE));

        // Content
        div.add(new Paragraph("The data on which this report is based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation.")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("The analysis done for this report has passed all internal quality controls.")
                .addStyle(BODY_TEXT_STYLE));

        div.add(new Paragraph("For feedback or complaints please contact ")
                .addStyle(BODY_TEXT_STYLE)
                .add(new Text(ReportResources.CONTACT_EMAIL_QA)
                        .addStyle(BODY_TEXT_BOLD_STYLE)));

        div.add(new Paragraph("For general questions, please contact us at ")
                .addStyle(BODY_TEXT_STYLE)
                .add(new Text(ReportResources.CONTACT_EMAIL_GENERAL)
                        .addStyle(BODY_TEXT_BOLD_STYLE)));

        return div;

    }


}
