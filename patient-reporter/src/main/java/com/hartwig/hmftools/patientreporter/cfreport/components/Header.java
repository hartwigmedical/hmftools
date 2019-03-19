package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;

public final class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg"; // "src/main/resources/pdf/hartwig_logo.jpg";
    private static PdfImageXObject hmfLogoObj = null;

    public static void addHeader(PdfPage page) {

        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        // Attempt to load image object
        if (hmfLogoObj == null) {
            hmfLogoObj = new PdfImageXObject(ReportResources.loadImageData(HMF_LOGO_PATH));
        }

        // Add image object to page
        if (hmfLogoObj != null) {
            canvas.addXObject(hmfLogoObj, 52, 772, 44, false);
        }

        // @TODO: add "Hartwig Medical OncoAct"


        canvas.release();

    }

}
