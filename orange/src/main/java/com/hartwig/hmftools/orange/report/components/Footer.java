package com.hartwig.hmftools.orange.report.components;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.TextAlignment;

import org.jetbrains.annotations.NotNull;

public class Footer
{
    private final List<FooterTemplate> footerTemplates = Lists.newArrayList();
    @NotNull
    private final ReportResources reportResources;
    private final boolean addDisclaimer;

    public Footer(@NotNull ReportResources reportResources, boolean addDisclaimer)
    {
        this.reportResources = reportResources;
        this.addDisclaimer = addDisclaimer;
    }

    public void renderFooter(@NotNull PdfPage page)
    {
        PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        int pageNumber = page.getDocument().getPageNumber(page);
        PdfFormXObject pageNumberTemplate = new PdfFormXObject(new Rectangle(0, 0, 465, 20));
        canvas.addXObject(pageNumberTemplate, 58, 20);
        footerTemplates.add(new FooterTemplate(pageNumber, pageNumberTemplate, addDisclaimer));

        canvas.release();
    }

    public void writeFooters(@NotNull PdfDocument document)
    {
        int totalPageCount = document.getNumberOfPages();
        for(FooterTemplate tpl : footerTemplates)
        {
            tpl.renderFooter(totalPageCount, document, reportResources);
        }
    }

    private static class FooterTemplate
    {

        private final int pageNumber;
        @NotNull
        private final PdfFormXObject template;
        private final boolean addDisclaimer;

        FooterTemplate(int pageNumber, @NotNull PdfFormXObject template, boolean addDisclaimer)
        {
            this.pageNumber = pageNumber;
            this.template = template;
            this.addDisclaimer = addDisclaimer;
        }

        void renderFooter(int totalPageCount, @NotNull PdfDocument document, @NotNull ReportResources reportResources)
        {
            String displayString = pageNumber + "/" + totalPageCount;

            Canvas canvas = new Canvas(template, document);
            Paragraph pageNumberParagraph = new Paragraph().add(displayString).addStyle(reportResources.pageNumberStyle());
            canvas.showTextAligned(pageNumberParagraph, 0, 0, TextAlignment.LEFT);

            if(addDisclaimer)
            {
                List<Text> disclaimerParts = List.of(
                        new Text("This report is for '"),
                        new Text("Research Use Only (RUO)").setFont(reportResources.fontBold()),
                        new Text("' and is "),
                        new Text("not").setUnderline(),
                        new Text(" suitable for diagnostic or clinical applications. "),
                        new Text("No rights").setUnderline(),
                        new Text(" can be derived from the contents of this report.")
                );
                Paragraph disclaimerParagraph = new Paragraph().setMaxWidth(430).addStyle(reportResources.disclaimerStyle());
                disclaimerParagraph.addAll(disclaimerParts);
                canvas.showTextAligned(disclaimerParagraph, 40, 0, TextAlignment.LEFT);
            }
        }
    }
}
