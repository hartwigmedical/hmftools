package com.hartwig.hmftools.orange.report.components;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public final class SidePanel {

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
        canvas.rectangle(pageSize.getWidth(), pageSize.getHeight(), -RECTANGLE_WIDTH, -RECTANGLE_HEIGHT);
        canvas.setFillColor(ReportResources.PALETTE_ORANGE);
        canvas.fill();

        int sideTextIndex = 0;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(pageSize, ++sideTextIndex, "Sample", report.sampleId()));
        cv.add(createSidePanelDiv(pageSize, ++sideTextIndex, "Platinum version", pipelineVersion(report)));

        canvas.release();
    }

    @NotNull
    private static String pipelineVersion(@NotNull OrangeReport report) {
        return report.platinumVersion() != null ? report.platinumVersion() : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private static Div createSidePanelDiv(@NotNull Rectangle pageSize, int index, @NotNull String label, @NotNull String value) {
        Div div = new Div();
        div.setKeepTogether(true);

        float yPos = (pageSize.getHeight() + 15) - index * ROW_SPACING;
        float xPos = pageSize.getWidth() - RECTANGLE_WIDTH + 15;

        div.add(new Paragraph(label.toUpperCase()).addStyle(ReportResources.sidePanelLabelStyle())
                .setFixedPosition(xPos, yPos, MAX_WIDTH));

        float valueFontSize = maxPointSizeForWidth(ReportResources.fontBold(), 11, 6, value, MAX_WIDTH);
        yPos -= VALUE_TEXT_Y_OFFSET;
        div.add(new Paragraph(value).addStyle(ReportResources.sidePanelValueStyle().setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(xPos, yPos, MAX_WIDTH)
                .setFixedLeading(valueFontSize));

        return div;
    }

    private static float maxPointSizeForWidth(@NotNull PdfFont font, float initialFontSize, float minFontSize, @NotNull String text,
            float maxWidth) {
        float fontIncrement = 0.1F;

        float fontSize = initialFontSize;
        float width = font.getWidth(text, initialFontSize);
        while (width > maxWidth && fontSize > minFontSize) {
            fontSize -= fontIncrement;
            width = font.getWidth(text, fontSize);
        }

        return fontSize;
    }
}
