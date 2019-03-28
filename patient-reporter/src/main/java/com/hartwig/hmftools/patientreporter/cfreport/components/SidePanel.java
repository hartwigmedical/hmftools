package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

import java.util.StringJoiner;

public final class SidePanel {

    private static final float ROW_SPACING = 42;
    private static final float CONTENT_X_START = 455;

    public static void renderSidePanel(PdfPage page, @NotNull AnalysedPatientReport patientReport, boolean fullHeight, boolean fullContent) {

        final SampleReport sampleReport = patientReport.sampleReport();

        // Draw background and markers
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.renderMarkerGrid(4, (fullHeight ? 20 : 2), CONTENT_X_START, 35, 820, -ROW_SPACING, .05f, .15f, canvas);

        // Add side panel content that is always on the side panel (full height or not)
        int sideTextIndex = 0;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(sideTextIndex++, "HMF sample id", patientReport.sampleReport().sampleId()));
        cv.add(createSidePanelDiv(sideTextIndex++, "Report date", ReportResources.REPORT_DATE));

        // Add side panel content that is only on the summary page
        if (fullHeight && fullContent) {

            if (!sampleReport.contactNames().isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Name requestor", sampleReport.contactNames()));
            }

            if (!sampleReport.contactEmails().isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Email requestor", sampleReport.contactEmails()));
            }

            // @TODO Hospital not in data?
//            if (!sampleReport.hospital().isEmpty()) {
//                cv.add(createSidePanelDiv(sideTextIndex++, "Hospital", sampleReport.hospital()));
//            }

            if (sampleReport.hospitalPatientId() != null && !sampleReport.hospitalPatientId().isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Hospital patient id", sampleReport.hospitalPatientId()));
            }

            // @TODO Patient gender not in data?
//            if (!sampleReport.patientGender().isEmpty()) {
//                cv.add(createSidePanelDiv(sideTextIndex++, "Hospital", sampleReport.patientGender()));
//            }

            // @TODO Patient birth date not in data?
//            if (sampleReport.patientBirthDate() != null) {
//                cv.add(createSidePanelDiv(sideTextIndex, "Hospital", ReportResources.formatDate(sampleReport.patientBirthDate())));
//            }

        }

        canvas.release();

    }

    /**
     * Draw background rectangle, either full height or only top of the page
     */
    private static void renderBackgroundRect(boolean fullHeight, @NotNull PdfCanvas canvas, @NotNull Rectangle pageSize) {

        final float RECTANGLE_WIDTH = 170;           // Width of the blue rectangle in pt
        final float RECTANGLE_HEIGHT_SHORT = 110;    // Height of the blue rectangle in pt when not full page height

        canvas.rectangle(pageSize.getWidth(), pageSize.getHeight(), -RECTANGLE_WIDTH, fullHeight ? -pageSize.getHeight() : -RECTANGLE_HEIGHT_SHORT);
        canvas.setFillColor(ReportResources.PALETTE_BLUE);
        canvas.fill();

    }

    @NotNull
    private static Div createSidePanelDiv(int index, @NotNull String label, @NotNull String value) {

        final float Y_START = 802;
        final float VALUE_TEXT_Y_OFFSET = 18;
        final float MAX_WIDTH = 120;

        Div div = new Div();
        div.setKeepTogether(true);

        // Add label
        float yPos = Y_START - index * ROW_SPACING;
        div.add(new Paragraph(label.toUpperCase())
                .addStyle(ReportResources.sidepanelLabelStyle())
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH));


        // Add value (auto resize the font if needed)
        final float valueFontSize = ReportResources.getMaxPointSizeForWidth(ReportResources.getFontBold(), 11, 6, value, MAX_WIDTH);
        yPos -= VALUE_TEXT_Y_OFFSET;
        div.add(new Paragraph(value)
                .addStyle(ReportResources.sidepanelValueStyle()
                        .setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH)
                .setFixedLeading(valueFontSize)
        );

        return div;

    }

 }
