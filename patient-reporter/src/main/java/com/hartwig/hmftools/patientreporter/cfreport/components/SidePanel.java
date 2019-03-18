package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportConfiguration;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;
import org.jetbrains.annotations.NotNull;

public class SidePanel {

    public enum PageMode {
        SummaryPage,
        ClosingPage,
        ContentPage
    };

    public static void draw(PageMode mode, @NotNull PdfWriter pdfWriter) {

        // Draw background and markers
        drawBackground(mode, pdfWriter);
        BaseMarker.drawSidepanelMarkers(mode == PageMode.SummaryPage || mode == PageMode.ClosingPage, pdfWriter);

        // @TODO Add sidepanel content from SampleReport

    }

    private static void drawBackground(PageMode mode, @NotNull PdfWriter writer) {

        PdfContentByte pb = writer.getDirectContent();
        pb.saveState();

        Rectangle pageSize = writer.getPageSize();

        pb.rectangle(pageSize.getWidth(), pageSize.getHeight(), -170, (mode == PageMode.SummaryPage || mode == PageMode.ClosingPage) ? -pageSize.getHeight() : -110);
        pb.setColorFill(ReportConfiguration.PALETTE_BLUE);
        pb.fill();
        pb.restoreState();

    }

}
