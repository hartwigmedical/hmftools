package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportConfiguration;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;
import org.jetbrains.annotations.NotNull;

public class SidePanel {

    public static void draw(boolean firstPage, @NotNull PdfWriter pdfWriter) {

        // Draw background and markers
        drawBackground(firstPage, pdfWriter);
        BaseMarker.drawSidepanelMarkers(firstPage, pdfWriter);

        // @TODO Add sidepanel content from SampleReport

    }

    private static void drawBackground(boolean fullheight, @NotNull PdfWriter writer) {

        PdfContentByte pb = writer.getDirectContent();
        pb.saveState();

        Rectangle pageSize = writer.getPageSize();

        pb.rectangle(pageSize.getWidth(), pageSize.getHeight(), -170, fullheight ? -pageSize.getHeight() : -110);
        pb.setColorFill(ReportConfiguration.PALETTE_BLUE);
        pb.fill();
        pb.restoreState();

    }

}
