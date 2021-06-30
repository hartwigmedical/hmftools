package com.hartwig.hmftools.orange.report.component;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
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
    private static final float RECTANGLE_WIDTH = 170;
    private static final float RECTANGLE_HEIGHT_SHORT = 110;

    private SidePanel() {
    }

    public static void renderSidePanel(@NotNull PdfPage page, @NotNull OrangeReport report) {
        PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(canvas, pageSize);

        int sideTextIndex = -1;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(++sideTextIndex, "Sample", report.sampleId()));
        cv.add(createSidePanelDiv(++sideTextIndex, "Pipeline", "v5.23"));

        canvas.release();
    }

    private static void renderBackgroundRect(@NotNull PdfCanvas canvas, @NotNull Rectangle pageSize) {
        canvas.rectangle(pageSize.getWidth(), pageSize.getHeight(), -RECTANGLE_WIDTH, -RECTANGLE_HEIGHT_SHORT);
        canvas.setFillColor(ReportResources.PALETTE_ORANGE);
        canvas.fill();
    }

    @NotNull
    private static Div createSidePanelDiv(int index, @NotNull String label, @NotNull String value) {
        float yStart = 802;
        float valueTextYOffset = 18;
        float maxWidth = 120;

        Div div = new Div();
        div.setKeepTogether(true);

        float yPos = yStart - index * ROW_SPACING;
        div.add(new Paragraph(label.toUpperCase()).addStyle(ReportResources.sidePanelLabelStyle())
                .setFixedPosition(CONTENT_X_START, yPos, maxWidth));

        float valueFontSize = ReportResources.maxPointSizeForWidth(ReportResources.fontBold(), 11, 6, value, maxWidth);
        yPos -= valueTextYOffset;
        div.add(new Paragraph(value).addStyle(ReportResources.sidePanelValueStyle().setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(CONTENT_X_START, yPos, maxWidth)
                .setFixedLeading(valueFontSize));

        return div;
    }
}
