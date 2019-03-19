package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.PageEventHandler;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.TextAlignment;

public final class Footer {

    // Total page count template
    private static final float PAGECOUNT_WIDTH = 20;
    private static final float PAGECOUNT_HEIGHT = 20;
    private static final float PAGECOUNT_X = 58;
    private static final float PAGECOUNT_Y = 20;
    private static final float PAGECOUNT_HSPACING = 0;
    private static final float PAGECOUNT_DESCENT = 0;
    private static final PdfFormXObject PAGECOUNT_PLACEHOLDER = new PdfFormXObject(new Rectangle(0, 0, PAGECOUNT_WIDTH, PAGECOUNT_HEIGHT));

    public static void addFooter(PdfPage page, PageEventHandler.PageMode pageMode) {

        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        // Add page number
        int pageNumber = page.getDocument().getPageNumber(page);
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());
        Paragraph p = new Paragraph().add(String.valueOf(pageNumber)).add("/");
        cv.showTextAligned(p, PAGECOUNT_X, PAGECOUNT_Y, TextAlignment.CENTER.RIGHT);
        canvas.addXObject(PAGECOUNT_PLACEHOLDER, PAGECOUNT_X + PAGECOUNT_HSPACING, PAGECOUNT_Y - PAGECOUNT_DESCENT);

        // Draw markers
        final boolean fullWidth = (pageMode != PageEventHandler.PageMode.SummaryPage && pageMode != PageEventHandler.PageMode.ClosingPage);
        BaseMarker.drawMarkerGrid(fullWidth ? 5 : 3,1,156, 87, 22, 0, .2f, 0, canvas);

        canvas.release();

    }

    public static void writeTotalPageCount(PdfDocument document) {
        Canvas canvas = new Canvas(PAGECOUNT_PLACEHOLDER, document);
        canvas.showTextAligned(String.valueOf(document.getNumberOfPages()), 0, PAGECOUNT_DESCENT, TextAlignment.LEFT);
    }

}
