package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.lowagie.text.Document;
import com.lowagie.text.pdf.PdfWriter;

public class Footer {

    public static void drawFooter(Document document, PdfWriter writer) {

        // @TODO Add page number

        // Footer line
        BaseMarker.drawFooterLine((document.getPageNumber() > 1), writer);

    }

}
