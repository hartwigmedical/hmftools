package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfPage;

public class PageEventHandler implements IEventHandler {

    public enum PageMode {
        SummaryPage,        // Full height sidepanel, full content
        ClosingPage,        // Full height sidepanel, just ID and report date
        ContentPage         // Short sidepanel with just ID and report date
    };

    private PageMode pageMode = PageMode.SummaryPage;

    public void setPageMode(PageMode pageMode) {
        this.pageMode = pageMode;
    }

    @Override
    public void handleEvent(Event event) {

        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType() == PdfDocumentEvent.START_PAGE) {

            final PdfPage page = documentEvent.getPage();

            // Add sidepanel and footer to new page
            Header.addHeader(page);
            SidePanel.addSidePanel(page, pageMode);
            Footer.addFooter(page, pageMode);

        }

    }

}
