package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.TextAlignment;
import org.jetbrains.annotations.NotNull;

public class Footer {

    // Total page count template
    private static final float PAGE_COUNT_WIDTH = 20;
    private static final float PAGE_COUNT_HEIGHT = 20;
    private static final float PAGE_COUNT_X = 58;
    private static final float PAGE_COUNT_Y = 20;
    private static final float PAGE_COUNT_HSPACING = .8f;
    private static final float PAGE_COUNT_DESCENT = 0;

    private PdfFormXObject pageCountPlaceholder = new PdfFormXObject(new Rectangle(0, 0, PAGE_COUNT_WIDTH, PAGE_COUNT_HEIGHT));

    public void renderFooter(PdfPage page, boolean fullWidth) {

        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        // Add current page number
        int pageNumber = page.getDocument().getPageNumber(page);
        Paragraph p = createPageNumberParagraph(String.format("%d/", pageNumber));
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());
        cv.showTextAligned(p, PAGE_COUNT_X, PAGE_COUNT_Y, TextAlignment.CENTER.RIGHT);

        // Add placeholder for total page count
        canvas.addXObject(pageCountPlaceholder, PAGE_COUNT_X + PAGE_COUNT_HSPACING, PAGE_COUNT_Y - PAGE_COUNT_DESCENT);

        // Draw markers
        BaseMarker.renderMarkerGrid(fullWidth ? 5 : 3,1,156, 87, 22, 0, .2f, 0, canvas);

        canvas.release();

    }

    public void writeTotalPageCount(@NotNull PdfDocument document) {
        Canvas canvas = new Canvas(pageCountPlaceholder, document);
        Paragraph p = createPageNumberParagraph(String.valueOf(document.getNumberOfPages()));
        canvas.showTextAligned(p, 0, PAGE_COUNT_DESCENT, TextAlignment.LEFT);
    }

    private static Paragraph createPageNumberParagraph(@NotNull String text) {
        return new Paragraph()
                .add(text)
                .addStyle(ReportResources.pageNumberStyle());
    }

}
