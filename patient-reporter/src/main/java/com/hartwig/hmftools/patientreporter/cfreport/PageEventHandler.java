package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfPage;

public class PageEventHandler implements IEventHandler {

    private boolean fullSidebar = true;
    private boolean fullSidebarContent = true;

    public void setSidebarType(boolean full, boolean fullContent) {
        fullSidebar = full;
        fullSidebarContent = fullSidebar ? fullContent : false;
    }

    @Override
    public void handleEvent(Event event) {

        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType() == PdfDocumentEvent.START_PAGE) {

            final PdfPage page = documentEvent.getPage();

            // Add sidepanel and footer to new page
            Header.addHeader(page);
            SidePanel.addSidePanel(page, fullSidebar, fullSidebarContent);
            Footer.addFooter(page, !fullSidebar);

        }

    }

}
