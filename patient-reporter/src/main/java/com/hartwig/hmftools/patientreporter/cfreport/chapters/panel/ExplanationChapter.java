package com.hartwig.hmftools.patientreporter.cfreport.chapters.panel;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.panel.PanelReport;
import com.itextpdf.io.IOException;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ExplanationChapter implements ReportChapter {

    @NotNull
    private final PanelReport report;

    public ExplanationChapter(@NotNull PanelReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @NotNull
    @Override
    public String name() {
        return "Result explanation";
    }

    @Override
    public void render(@NotNull Document reportDocument) throws IOException {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createExplanationDiv()));
        reportDocument.add(table);
    }

    @NotNull
    private static Div createExplanationDiv() {
        Div div = new Div();

        div.add(new Paragraph("Details on the Result File").addStyle(ReportResources.smallBodyHeadingStyle()));

        div.add(createContentParagraph("The variant calling of the sequencing data is based on reference genome version GRCh37."));
        div.add(createContentParagraph("Transcript list can be found on."));
        div.add(createContentParagraph("Short description of the headers present in the VCF"));
        return div;
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }
}