package com.hartwig.hmftools.orange.report;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    @NotNull
    private final String outputDir;

    public ReportWriter(@NotNull final String outputDir) {
        this.outputDir = outputDir;
    }

    public void write(@NotNull OrangeReport report, boolean writeToFile) throws IOException {
        LOGGER.info("Writing {}", report);
        String outputFilePath = outputDir + File.separator + report.sampleId() + ".orange.pdf";

        Document doc = initializeReport(outputFilePath, writeToFile);

        doc.getPdfDocument().addNewPage();
        renderTitle(doc.getPdfDocument().getFirstPage());
        Image circosImage = new Image(ImageDataFactory.create(report.plots().purpleCircosPlot()));
        doc.add(circosImage);

        doc.close();

        if (writeToFile) {
            LOGGER.info("Created ORANGE report at {}", outputFilePath);
        } else {
            LOGGER.info("Successfully generated in-memory ORANGE report");
        }
    }

    private void renderTitle(@NotNull PdfPage page) {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        cv.add(new Paragraph().add(new Text("ORANGE Report").setFontSize(11).setFixedPosition(230, 791, 300)));

        PdfFormXObject chapterTitleTemplate = new PdfFormXObject(new Rectangle(0, 0, 500, 30));
        pdfCanvas.addXObject(chapterTitleTemplate, 29, 721);

        pdfCanvas.release();
    }

    @NotNull
    private static Document initializeReport(@NotNull String outputFilePath, boolean writeToFile) throws IOException {
        PdfWriter writer;
        if (writeToFile) {
            writer = new PdfWriter(outputFilePath);
        } else {
            // Write output to output stream where it is effectively ignored.
            writer = new PdfWriter(new ByteArrayOutputStream());
        }

        PdfDocument pdf = new PdfDocument(writer);
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle("ORANGE Report");
        pdf.getDocumentInfo().setAuthor("Platinum");

        Document document = new Document(pdf);
        document.setMargins(62, 29, 62, 55.5F);

        return document;
    }
}
