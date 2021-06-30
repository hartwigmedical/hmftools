package com.hartwig.hmftools.orange.report.component;

import java.net.MalformedURLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Header {

    private static final Logger LOGGER = LogManager.getLogger(Header.class);

    @Nullable
    private final PdfImageXObject orangeImageObj;
    private final List<ChapterPageCounter> chapterPageCounters = Lists.newArrayList();

    public Header(@NotNull String orangeImagePath) {
        ImageData orangeImage = null;
        try {
            orangeImage = ImageDataFactory.create(orangeImagePath);
        } catch (MalformedURLException e) {
            LOGGER.warn("Could not orange image from {}", orangeImagePath);
        }
        orangeImageObj = orangeImage != null ? new PdfImageXObject(orangeImage) : null;
    }

    public void renderHeader(@NotNull String chapterTitle, boolean firstPageOfChapter, @NotNull PdfPage page) {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        if (orangeImageObj != null) {
            pdfCanvas.addXObject(orangeImageObj, 52, 772, 60, false);
        }

        cv.add(new Paragraph().add(new Text("O").setFont(ReportResources.fontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_ORANGE_1))
                .add(new Text("R").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_2))
                .add(new Text("A").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_3))
                .add(new Text("N").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_4))
                .add(new Text("G").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_5))
                .add(new Text("E").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_6))
                .add(new Text(" Report").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLACK))
                .setFixedPosition(200, 800, 300));

        if (firstPageOfChapter) {
            chapterPageCounters.add(new ChapterPageCounter(chapterTitle));
        }

        PdfFormXObject chapterTitleTemplate = new PdfFormXObject(new Rectangle(0, 0, 500, 30));
        pdfCanvas.addXObject(chapterTitleTemplate, ReportResources.PAGE_MARGIN_LEFT, 721);
        chapterPageCounters.get(chapterPageCounters.size() - 1).addPage(chapterTitleTemplate);

        pdfCanvas.release();
    }

    public void writeChapterTitles(@NotNull PdfDocument document) {
        for (ChapterPageCounter cpc : chapterPageCounters) {
            cpc.renderChapterTitles(document);
        }
    }

    private static class ChapterPageCounter {

        @NotNull
        private final String chapterTitle;
        private final List<PdfFormXObject> templates = Lists.newArrayList();

        ChapterPageCounter(@NotNull String chapterTitle) {
            this.chapterTitle = chapterTitle;
        }

        void addPage(@NotNull PdfFormXObject chapterTitleTemplate) {
            templates.add(chapterTitleTemplate);
        }

        void renderChapterTitles(@NotNull PdfDocument document) {
            int totalChapterPages = templates.size();

            for (int i = 0; i < templates.size(); i++) {
                PdfFormXObject tpl = templates.get(i);

                String text = chapterTitle + (totalChapterPages > 1 ? " (" + (i + 1) + "/" + totalChapterPages + ")" : "");

                Canvas canvas = new Canvas(tpl, document);
                canvas.add(new Paragraph(text).addStyle(ReportResources.chapterTitleStyle()));
            }
        }
    }
}
