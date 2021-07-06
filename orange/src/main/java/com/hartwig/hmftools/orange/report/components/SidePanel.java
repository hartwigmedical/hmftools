package com.hartwig.hmftools.orange.report.components;

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

    private static final float CONTENT_X_START = 455;
    private static final float CONTENT_Y_START = 820;
    private static final float ROW_SPACING = 35;
    private static final float VALUE_TEXT_Y_OFFSET = 18;
    private static final float MAX_WIDTH = 120;

    private static final float RECTANGLE_WIDTH = 170;
    private static final float RECTANGLE_HEIGHT = 84;

    private SidePanel() {
    }

    public static void renderSidePanel(@NotNull PdfPage page, @NotNull OrangeReport report) {
        PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(canvas, pageSize);

        int sideTextIndex = -1;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(++sideTextIndex, "Sample", report.sampleId()));
        cv.add(createSidePanelDiv(++sideTextIndex, "Pipeline", report.pipelineVersion()));

        canvas.release();
    }

    private static void renderBackgroundRect(@NotNull PdfCanvas canvas, @NotNull Rectangle pageSize) {
        canvas.rectangle(pageSize.getWidth(), pageSize.getHeight(), -RECTANGLE_WIDTH, -RECTANGLE_HEIGHT);
        canvas.setFillColor(ReportResources.PALETTE_ORANGE);
        canvas.fill();
    }

    @NotNull
    private static Div createSidePanelDiv(int index, @NotNull String label, @NotNull String value) {
        Div div = new Div();
        div.setKeepTogether(true);

        float yPos = CONTENT_Y_START - index * ROW_SPACING;
        div.add(new Paragraph(label.toUpperCase()).addStyle(ReportResources.sidePanelLabelStyle())
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH));

        float valueFontSize = ReportResources.maxPointSizeForWidth(ReportResources.fontBold(), 11, 6, value, MAX_WIDTH);
        yPos -= VALUE_TEXT_Y_OFFSET;
        div.add(new Paragraph(value).addStyle(ReportResources.sidePanelValueStyle().setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH)
                .setFixedLeading(valueFontSize));

        return div;
    }
}
