package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;

import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.lowagie.text.Document;
import com.lowagie.text.pdf.PdfPageEventHelper;
import com.lowagie.text.pdf.PdfWriter;

public class PageEvent extends PdfPageEventHelper {

    private SidePanel.PageMode sidePanelPageMode = SidePanel.PageMode.SummaryPage;

    public void setSidePanelPageMode(SidePanel.PageMode pageMode) {
        sidePanelPageMode = pageMode;
    }


    @Override
    public void onStartPage(PdfWriter writer, Document document) {

        //@TODO Add header

        // Add side panel
        SidePanel.draw(sidePanelPageMode, writer);

        // Add footer
        Footer.drawFooter(document, writer);

    }

    @Override
    public void onCloseDocument(PdfWriter writer, Document document) { }

}

