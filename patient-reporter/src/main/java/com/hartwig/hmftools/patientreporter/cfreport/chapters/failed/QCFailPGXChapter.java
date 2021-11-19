package com.hartwig.hmftools.patientreporter.cfreport.chapters.failed;

import java.util.List;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.Pharmacogenetics;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class QCFailPGXChapter implements ReportChapter {

    @NotNull
    private final QCFailReport failReport;

    public QCFailPGXChapter(@NotNull QCFailReport failReport) {
        this.failReport = failReport;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    public String name() {
        return "pharmagenetics";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        reportDocument.add(createPeachGenotypesTable(failReport.peachGenotypes(), failReport.sampleReport().reportPharmogenetics()));

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        table.setWidth(contentWidth());

        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported pharmacogenetics")));
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDivWithLinkThree("See the directory 'Patient Reporting' in ",
                        "https://resources.hartwigmedicalfoundation.nl ",
                        "for details on the panel and for more links to advice on treatment adjustments.",
                        "https://resources.hartwigmedicalfoundation.nl"))
                .add(createContentDiv(new String[] {
                        "The called haplotypes for a gene are the simplest combination of haplotypes that perfectly explains all of the "
                                + "observed variants for that gene. If no combination of haplotypes in the panel can perfectly explain the "
                                + "observed variants, then 'Unresolved Haplotype' is called.",
                        "Wild type is assumed when no variants are observed." })));
        table.addCell(TableUtil.createLayoutCell());
        reportDocument.add(table);

    }

    @NotNull
    private static Table createPeachGenotypesTable(@NotNull List<PeachGenotype> peachGenotypes, boolean reportPeach) {
        String title = "Pharmacogenetics";

        if (reportPeach) {
            if (peachGenotypes.isEmpty()) {
                return TableUtil.createNoneReportTable(title);
            } else  {
                Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 60, 60, 100, 60 },
                        new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Genotype"),
                                TableUtil.createHeaderCell("Function"), TableUtil.createHeaderCell("Linked drugs"),
                                TableUtil.createHeaderCell("Source").setTextAlignment(TextAlignment.CENTER) });

                for (PeachGenotype peachGenotype : Pharmacogenetics.sort(peachGenotypes)) {
                    contentTable.addCell(TableUtil.createContentCell(peachGenotype.gene()));
                    contentTable.addCell(TableUtil.createContentCell(peachGenotype.haplotype()));
                    contentTable.addCell(TableUtil.createContentCell(peachGenotype.function()));
                    contentTable.addCell(TableUtil.createContentCell(peachGenotype.linkedDrugs()));
                    contentTable.addCell(TableUtil.createContentCell(new Paragraph(Pharmacogenetics.sourceName(peachGenotype.urlPrescriptionInfo())).addStyle(
                                    ReportResources.dataHighlightLinksStyle()))
                            .setAction(PdfAction.createURI(Pharmacogenetics.url(peachGenotype.urlPrescriptionInfo())))
                            .setTextAlignment(TextAlignment.CENTER));
                }
                return TableUtil.createWrappingReportTable(title, contentTable);
            }
        } else {
            return TableUtil.createNAReportTable(title);
        }
    }

    @NotNull
    private static Div createContentDiv(@NotNull String[] contentParagraphs) {
        Div div = new Div();
        for (String s : contentParagraphs) {
            div.add(new Paragraph(s).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING));
        }
        return div;
    }

    @NotNull
    private static Paragraph createSectionTitle(@NotNull String sectionTitle) {
        return new Paragraph(sectionTitle).addStyle(ReportResources.smallBodyHeadingStyle());
    }

    @NotNull
    private static Paragraph createParaGraphWithLinkThree(@NotNull String string1, @NotNull String string2, @NotNull String string3,
            @NotNull String link) {
        return new Paragraph(string1).addStyle(ReportResources.subTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string2).addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(link)))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string3).addStyle(ReportResources.subTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Div createContentDivWithLinkThree(@NotNull String string1, @NotNull String string2, @NotNull String string3,
            @NotNull String link) {
        Div div = new Div();

        div.add(createParaGraphWithLinkThree(string1, string2, string3, link));
        return div;
    }
}