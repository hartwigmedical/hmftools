package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

public final class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg";

    private PdfImageXObject hmfLogoObj;

    private ArrayList<ChapterPageCounter> chapterPageCounters = new ArrayList<>();

    public Header() {

        // Attempt to load image object
        ImageData imgData = ReportResources.loadImageData(HMF_LOGO_PATH);
        if (imgData != null) {
            hmfLogoObj = new PdfImageXObject(imgData);
        }

    }

    public void renderHeader(@NotNull String chapterTitle, boolean firstPageOfChapter, @NotNull PdfPage page) {

        final PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        // Add HMF logo image object to page
        if (hmfLogoObj != null) {
            pdfCanvas.addXObject(hmfLogoObj, 52, 772, 44, false);
        }

        // Add "Hartwig Medical OncoAct"
        cv.add(new Paragraph()
                .add(new Text("Hartwig Medical")
                        .setFont(ReportResources.getFontBold())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text(" Onco")
                        .setFont(ReportResources.getFontRegular())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text("Act")
                        .setFont(ReportResources.getFontRegular())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_RED))
                .setFixedPosition(230, 791, 300));


        if (firstPageOfChapter) {
            chapterPageCounters.add(new ChapterPageCounter(chapterTitle));
        }

        // Create chapter title template
        PdfFormXObject chapterTitleTemplate = new PdfFormXObject(new Rectangle(0, 0, 500, 20));
        pdfCanvas.addXObject(chapterTitleTemplate,ReportResources.PAGE_MARGIN_LEFT, 721);
        chapterPageCounters.get(chapterPageCounters.size() - 1).addPage(chapterTitleTemplate);

        pdfCanvas.release();

    }

    public void writeChapterTitles(@NotNull PdfDocument document) {
        for (ChapterPageCounter cpc: chapterPageCounters) {
            cpc. renderChapterTitles(document);
        }

    }

    public static class ChapterPageCounter {

        private String chapterTitle;
        private ArrayList<PdfFormXObject> templates = new ArrayList<>();

        public ChapterPageCounter(String chapterTitle) {
            this.chapterTitle = chapterTitle;
        }

        public void addPage(@NotNull PdfFormXObject chapterTitleTemplate) {
            templates.add(chapterTitleTemplate);
        }

        public void renderChapterTitles(@NotNull PdfDocument document) {

            int totalChapterPages = templates.size();

            for (int i = 0; i < templates.size(); i++) {
                PdfFormXObject tpl = templates.get(i);

                String text = chapterTitle
                        + (totalChapterPages > 1 ? " (" + (i+1) + "/" + totalChapterPages + ")" : "");

                Canvas canvas = new Canvas(tpl, document);
                canvas.add(new Paragraph(text).addStyle(ReportResources.chapterTitleStyle()));
            }
        }

    }

}
