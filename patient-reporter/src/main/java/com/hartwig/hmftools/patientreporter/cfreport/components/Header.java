package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;

import java.io.InputStream;
import java.net.MalformedURLException;

public class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg"; // "src/main/resources/pdf/hartwig_logo.jpg";
    private static PdfImageXObject hmfLogoObj = null;

    public static void addHeader(PdfPage page) {

        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        // Add image object to page
        if (hmfLogoObj != null) {
            canvas.addXObject(hmfLogoObj, 52, 772, 44, false);
        }

        // @TODO: add "Hartwig Medical OncoAct"


        canvas.release();

    }

    public static void loadHmfLogo(ClassLoader classLoader) throws java.io.IOException {

        // Load image data
        byte[] data = null;
        InputStream is = classLoader.getResourceAsStream(HMF_LOGO_PATH);
        data = new byte[is.available()];
        is.read(data);

        // Create Pdf XObject
        if (data != null) {
            ImageData imgData = ImageDataFactory.create(data, true);
            Header.hmfLogoObj = new PdfImageXObject(imgData);
        }

    }

}
