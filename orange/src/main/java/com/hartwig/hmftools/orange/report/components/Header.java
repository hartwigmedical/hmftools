package com.hartwig.hmftools.orange.report.components;

import java.net.URL;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

public class Header
{
    private final PdfImageXObject mOrangeCircosObject;
    private final ReportResources mReportResources;
    private final boolean mAddDisclaimer;

    public Header(final URL orangeCircosPath, final ReportResources reportResources, boolean addDisclaimer)
    {
        mOrangeCircosObject = new PdfImageXObject(ImageDataFactory.create(orangeCircosPath));
        mReportResources = reportResources;
        mAddDisclaimer = addDisclaimer;
    }

    public void renderHeader(final PdfPage page)
    {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas canvas = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        pdfCanvas.addXObject(mOrangeCircosObject, 50, page.getPageSize().getHeight() - 70, 60, false);

        Paragraph title = new Paragraph().add(new Text("O").setFont(mReportResources.fontBold())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_ORANGE_1))
                .add(new Text("R").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_2))
                .add(new Text("A").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_3))
                .add(new Text("N").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_4))
                .add(new Text("G").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_5))
                .add(new Text("E").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_6))
                .add(new Text(" Report").setFont(mReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLACK));

        if(mAddDisclaimer)
        {
            title = title.add(new Text(" (Research Use Only)").setFont(mReportResources.fontBold())
                    .setFontSize(9)
                    .setFontColor(ReportResources.PALETTE_BLACK));
        }

        float left = mAddDisclaimer ? 150 : 180;
        float width = mAddDisclaimer ? 370 : 300;

        canvas.add(title.setFixedPosition(left, page.getPageSize().getHeight() - 40, width));

        pdfCanvas.release();
    }
}
