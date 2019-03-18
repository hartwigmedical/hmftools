package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.PageEventHandler;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;

public class Footer {

    public static void addFooter(PdfPage page, PageEventHandler.PageMode pageMode) {

        // @TODO Add page number

        // Draw markers
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final boolean fullWidth = (pageMode != PageEventHandler.PageMode.SummaryPage && pageMode != PageEventHandler.PageMode.ClosingPage);
        BaseMarker.drawMarkerGrid(fullWidth ? 5 : 3,1,156, 87, 22, 0, .2f, 0, canvas);
        canvas.release();

    }

}
