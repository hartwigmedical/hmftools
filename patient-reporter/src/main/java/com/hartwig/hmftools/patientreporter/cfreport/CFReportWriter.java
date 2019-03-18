package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;
import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.AreaBreakType;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import sun.misc.IOUtils;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.nio.file.Files;

public class CFReportWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    @Override
    public void writeAnalysedPatientReport(@NotNull final AnalysedPatientReport report, @NotNull final String outputFilePath) {
        // TODO!
    }

    @Override
    public void writeQCFailReport(@NotNull final QCFailReport report, @NotNull final String outputFilePath) {
        // Does not need to be implemented
    }

    public void createReport(String outputFilePath) {

        try {

            final Document report = initializeReport(outputFilePath);
            final PdfDocument document = report.getPdfDocument();

            // Setup page event handling (used for automatic generation of header, sidepanel and footer)
            PageEventHandler pageEvents = new PageEventHandler();
            pageEvents.setPageMode(PageEventHandler.PageMode.SummaryPage);
            document.addEventHandler(PdfDocumentEvent.START_PAGE, pageEvents);

            // Load resources (images and fonts)
            Header.loadHmfLogo(this.getClass().getClassLoader());

            // Add summary page
            report.add(new Paragraph("Summary"));

            // Add intermediate pages
            pageEvents.setPageMode(PageEventHandler.PageMode.ContentPage);
            report.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            report.add(new Paragraph("Therapy details"));

            // Add closing page
            pageEvents.setPageMode(PageEventHandler.PageMode.ClosingPage);
            report.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            report.add(new Paragraph("Sample Details & Disclaimer"));

            Footer.writeTotalPageCount(document);
            report.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Temporary main method for quick running while report is in development
     *
     * @TODO: Remove this
     * @deprecated
     */
    public static void main(String[] args) {

        CFReportWriter writer = new CFReportWriter();
        writer.createReport("/Users/Wilco/hmf/tmp/temp_" + String.valueOf(System.currentTimeMillis()) + ".pdf");

    }

    /**
     * Initialize report document with relevant metadata
     * @return
     */
    @NotNull
    private static Document initializeReport(@NotNull final String outputFilePath) throws IOException {

        // Prevent overwriting existing files
        if (Files.exists(new File(outputFilePath).toPath())) {
            throw new IOException("Could not write " + outputFilePath + " as it already exists.");
        }

        // Create PDF with metadata
        final PdfDocument pdf = new PdfDocument(new PdfWriter(outputFilePath));
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle(ReportConfiguration.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportConfiguration.METADATA_AUTHOR);

        // Create document
        final Document document = new Document(pdf);
        document.setMargins(ReportConfiguration.PAGE_MARGIN_LEFT,
                ReportConfiguration.PAGE_MARGIN_RIGHT,
                ReportConfiguration.PAGE_MARGIN_TOP,
                ReportConfiguration.PAGE_MARGIN_BOTTOM);

        return document;
    }

}
