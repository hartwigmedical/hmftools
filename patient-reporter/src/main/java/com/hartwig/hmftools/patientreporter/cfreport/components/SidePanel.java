package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

public final class SidePanel {

    private static final float ROW_SPACING = 42;
    private static final float CONTENT_X_START = 455;

    public static void renderSidePanel(PdfPage page, boolean fullHeight, boolean fullContent) {

        // Draw background and markers
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.renderMarkerGrid(4, (fullHeight ? 20 : 2), CONTENT_X_START, 35, 820, -ROW_SPACING, .05f, .15f, canvas);

        // Add sidepanel content that is always on the sidepoanel (full height or not)
        int sideTextIndex = 0;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(sideTextIndex++, "HMF sample id", "CORE-01-99-1234T"));
        cv.add(createSidePanelDiv(sideTextIndex++, "Report date", ReportResources.REPORT_DATE));

        // Add sidepanel content that is only on the summary page
        if (fullHeight && fullContent) {
            cv.add(createSidePanelDiv(sideTextIndex++, "Name requestor", "Dr. Nola pluijmen"));
            cv.add(createSidePanelDiv(sideTextIndex++, "Email requestor", "NolaPluijmen415@undefinedmedicalcenter.com"));
            cv.add(createSidePanelDiv(sideTextIndex++, "Hospital", "OLVG Oost"));
            cv.add(createSidePanelDiv(sideTextIndex++, "Hospital patiend id", "839493929"));
            cv.add(createSidePanelDiv(sideTextIndex++, "Gender", "Female"));
            cv.add(createSidePanelDiv(sideTextIndex++, "Birth date", "11-Nov-1954"));
        }

        canvas.release();

    }

    /**
     * Draw background rectangle, either full height or only top
     *
     * @param fullHeight
     * @param canvas
     * @param pageSize
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
        final float VALUE_TEXT_Y_OFFSET = 16;
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
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH));

        return div;

    }



 }
