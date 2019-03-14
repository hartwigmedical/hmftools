package com.hartwig.hmftools.patientreporter.cfreport;

import com.lowagie.text.Document;
import com.lowagie.text.pdf.PdfPageEventHelper;
import com.lowagie.text.pdf.PdfWriter;

public class PageEvent extends PdfPageEventHelper {

    public void onStartPage(PdfWriter writer, Document document) {
        //@TODO Add header
    }

    public void onEndPage(PdfWriter writer, Document document) {
        //@TODO Add footer
    }

}
