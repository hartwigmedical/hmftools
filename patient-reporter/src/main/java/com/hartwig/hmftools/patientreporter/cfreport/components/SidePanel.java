package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

import java.time.LocalDate;

public final class SidePanel {

    private static final float ROW_SPACING = 42;
    private static final float CONTENT_X_START = 455;

    public static void renderSidePanel(PdfPage page, @NotNull final SampleReport sampleReport, boolean fullHeight, boolean fullContent) {

        // Draw background and markers
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.renderMarkerGrid(4, (fullHeight ? 20 : 2), CONTENT_X_START, 35, 820, -ROW_SPACING, .05f, .15f, canvas);

        // Add side panel content that is always on the side panel (full height or not)
        int sideTextIndex = 0;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(sideTextIndex++, "HMF sample id", sampleReport.sampleId()));
        cv.add(createSidePanelDiv(sideTextIndex++, "Report date", ReportResources.REPORT_DATE));

        // Add side panel content that is only on the summary page
        if (fullHeight && fullContent) {

            final String contactNames = "Dr. Nola Pluijmen"; // @TODO Replace with sampleReport.contactNames() which can be null or empty string
            if (contactNames != null && !contactNames.isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Name requestor", contactNames));
            }

            final String contactEmails = "NolaPluijmen415@gmail.com"; // @TODO Replace with sampleReport.contactEmails() which can be null or empty string
            if (contactEmails != null && !contactEmails.isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Email requestor", contactEmails));
            }

            final String hospitalName = "OLVG Oost"; // @TODO Replace with sampleReport.hospital() which can be null or empty string
            if (hospitalName != null && !hospitalName.isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Hospital", hospitalName));
            }

            final String hospitalPatientId = sampleReport.hospitalPatientId();
            if (hospitalPatientId != null && !hospitalPatientId.isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Hospital patient id", hospitalPatientId));
            }

            final String patientGender = "Female"; // @TODO Replace with sampleReport.patientGender() which can be null or empty string
            if (patientGender != null && !patientGender.isEmpty()) {
                cv.add(createSidePanelDiv(sideTextIndex++, "Gender", patientGender));
            }

            final LocalDate patientBirthDate = LocalDate.of(1973, 10, 4); // @TODO Replace with sampleReport.patientBirthDate() which can be null
            if (patientBirthDate != null) {
                cv.add(createSidePanelDiv(sideTextIndex, "Birth date", DataUtil.formatDate(patientBirthDate)));
            }

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
