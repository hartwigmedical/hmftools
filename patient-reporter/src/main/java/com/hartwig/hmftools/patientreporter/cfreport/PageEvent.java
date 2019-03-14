package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;

import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.lowagie.text.Document;
import com.lowagie.text.pdf.PdfPageEventHelper;
import com.lowagie.text.pdf.PdfWriter;

public class PageEvent extends PdfPageEventHelper {

    public void onStartPage(PdfWriter writer, Document document) {

        //@TODO Add header

        // Add side panel
        SidePanel.draw((document.getPageNumber() == 1), writer);

        // Add footer
        Footer.drawFooter(document.getPageNumber(), writer);

    }

}
