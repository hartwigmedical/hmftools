package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;
import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableCell;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHeaderCell;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.DottedBorder;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.AreaBreakType;
import com.itextpdf.layout.property.UnitValue;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.awt.*;
import java.io.File;
import java.io.IOException;
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

            // Add summary page
//            report.add(new Paragraph("Summary"));

            // Add intermediate pages
            pageEvents.setPageMode(PageEventHandler.PageMode.ContentPage);
//            report.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            report.add(new Paragraph("Therapy details"));

            Table table = new Table(UnitValue.createPercentArray(new float[] {18, 12, 61, 8, 15, 6}));
            table.setWidth(509);

            // Add header
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("DRIVERS")));
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("MATCH")));
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("TREATMENTS")));
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("LEVEL OF EVIDENCE")));
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("RESPONSE")));
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph("SOURCE")));

            // Add content
            for (int i = 0; i < 40; i++) {
                table.addCell(new TableCell().add(new Paragraph("BRAF p.Val600Glu")));
                table.addCell(new TableCell().add(new Paragraph("Specific")));
                table.addCell(new TableCell().add(new Paragraph("Binimetinib + Encorafenib")));
                table.addCell(new TableCell().add(new Paragraph("A")));
                table.addCell(new TableCell().add(new Paragraph("Responsive")));
                table.addCell(new TableCell().add(new Paragraph("OncoKB")));
            }

            report.add(table);

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
        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);

        // Create document
        final Document document = new Document(pdf);
        document.setMargins(ReportResources.PAGE_MARGIN_LEFT,
                ReportResources.PAGE_MARGIN_RIGHT,
                ReportResources.PAGE_MARGIN_TOP,
                ReportResources.PAGE_MARGIN_BOTTOM);

        return document;
    }

}
