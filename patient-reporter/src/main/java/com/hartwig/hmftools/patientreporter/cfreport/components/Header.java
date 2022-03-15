package com.hartwig.hmftools.patientreporter.cfreport.components;

import java.net.MalformedURLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
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
    private final PdfImageXObject companyLogoObj;
    private final List<ChapterPageCounter> chapterPageCounters = Lists.newArrayList();

    public Header(@NotNull String logoCompanyPath) {
        ImageData companyLogoImage = null;
        try {
            companyLogoImage = ImageDataFactory.create(logoCompanyPath);
        } catch (MalformedURLException e) {
            LOGGER.warn("Could not load company logo image from {}", logoCompanyPath);
        }
        companyLogoObj = companyLogoImage != null ? new PdfImageXObject(companyLogoImage) : null;
    }

    public void renderHeader(@NotNull String chapterTitle, @NotNull String pdfTitle, boolean firstPageOfChapter, @NotNull PdfPage page,
            boolean isWgsReport) {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        if (companyLogoObj != null) {
            pdfCanvas.addXObject(companyLogoObj, 52, 772, 44, false);
        }

        cv.add(new Paragraph().add(new Text("Hartwig Medical").setFont(ReportResources.fontBold())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text(" Onco").setFont(ReportResources.fontRegular()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text(isWgsReport ? "Act" : "Panel").setFont(ReportResources.fontRegular())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_RED))
                .setFixedPosition(230, 791, 300));

        if (firstPageOfChapter) {
            chapterPageCounters.add(new ChapterPageCounter(chapterTitle, pdfTitle));
        }

        PdfFormXObject chapterTitleTemplate;
        if (pdfTitle.isEmpty()) {
            chapterTitleTemplate = new PdfFormXObject(new Rectangle(0, 0, 500, 30));
        } else {
            chapterTitleTemplate = new PdfFormXObject(new Rectangle(0, 0, 500, 60));
        }

        pdfCanvas.addXObject(chapterTitleTemplate, ReportResources.PAGE_MARGIN_LEFT, 710);
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
        private final String pdfTitle;
        @NotNull
        private final String chapterTitle;
        private final List<PdfFormXObject> templates = Lists.newArrayList();

        ChapterPageCounter(@NotNull String chapterTitle, @NotNull String pdfTitle) {
            this.chapterTitle = chapterTitle;
            this.pdfTitle = pdfTitle;
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

                if (!pdfTitle.isEmpty()) {
                    canvas.add(new Paragraph(pdfTitle).addStyle(ReportResources.chapterTitleStyle()));
                }

                canvas.add(new Paragraph(text).addStyle(ReportResources.chapterTitleStyle()));
            }
        }
    }
}
