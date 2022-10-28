package com.hartwig.hmftools.patientreporter.cfreport.chapters.failed;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
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
        return "Pharmacogenetics";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        reportDocument.add(createPharmacogeneticsGenotypesTable(failReport.pharmacogeneticsGenotypes(), failReport.sampleReport().reportPharmogenetics()));

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        table.setWidth(contentWidth());

        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported pharmacogenetics")));
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDivWithLinkThree(
                        "The details on the pharmacogenetics haplotypes and advice on related treatment adjustments can be downloaded from ",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip",
                        ".",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip"))
                .add(createContentDiv(new String[] {
                        "The called haplotypes for a gene are the simplest combination of haplotypes that perfectly explains all of the "
                                + "observed variants for that gene. If no combination of haplotypes in the panel can perfectly explain the "
                                + "observed variants, then 'Unresolved Haplotype' is called.",
                        "Wild type is assumed when no variants are observed." })));
        table.addCell(TableUtil.createLayoutCell());
        reportDocument.add(table);

    }

    @NotNull
    private static Table createPharmacogeneticsGenotypesTable(@NotNull Map<String, List<PeachGenotype>> pharmacogeneticsGenotypes, boolean reportPharmacogenetics) {
        String title = "Pharmacogenetics";

        if (reportPharmacogenetics) {
            if (pharmacogeneticsGenotypes.isEmpty()) {
                return TableUtil.createNoneReportTable(title, null, TableUtil.TABLE_BOTTOM_MARGIN, ReportResources.CONTENT_WIDTH_WIDE);
            } else {
                Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 60, 60, 100, 60 },
                        new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Genotype"),
                                TableUtil.createHeaderCell("Function"), TableUtil.createHeaderCell("Linked drugs"),
                                TableUtil.createHeaderCell("Source").setTextAlignment(TextAlignment.CENTER) },
                        ReportResources.CONTENT_WIDTH_WIDE);

                Set<String> sortedPharmacogenetics = Sets.newTreeSet(pharmacogeneticsGenotypes.keySet().stream().collect(Collectors.toSet()));
                for (String sortPharmacogenetics : sortedPharmacogenetics) {
                    List<PeachGenotype> pharmacogeneticsGenotypeList = pharmacogeneticsGenotypes.get(sortPharmacogenetics);
                    contentTable.addCell(TableUtil.createContentCell(sortPharmacogenetics));

                    Table tableGenotype = new Table(new float[] { 1 });
                    Table tableFunction = new Table(new float[] { 1 });
                    Table tableLinkedDrugs = new Table(new float[] { 1 });
                    Table tableSource = new Table(new float[] { 1 });

                    for (PeachGenotype peachGenotype : pharmacogeneticsGenotypeList) {
                        tableGenotype.addCell(TableUtil.createTransparentCell(peachGenotype.haplotype()));
                        tableFunction.addCell(TableUtil.createTransparentCell(peachGenotype.function()));
                        tableLinkedDrugs.addCell(TableUtil.createTransparentCell(peachGenotype.linkedDrugs()));
                        tableSource.addCell(TableUtil.createTransparentCell(new Paragraph(Pharmacogenetics.sourceName(peachGenotype.urlPrescriptionInfo())).addStyle(
                                        ReportResources.dataHighlightLinksStyle()))
                                .setAction(PdfAction.createURI(Pharmacogenetics.url(peachGenotype.urlPrescriptionInfo()))));
                    }

                    contentTable.addCell(TableUtil.createContentCell(tableGenotype));
                    contentTable.addCell(TableUtil.createContentCell(tableFunction));
                    contentTable.addCell(TableUtil.createContentCell(tableLinkedDrugs));
                    contentTable.addCell(TableUtil.createContentCell(tableSource));
                }
                return TableUtil.createWrappingReportTable(title, null, contentTable, TableUtil.TABLE_BOTTOM_MARGIN);
            }
        } else {
            String noConsent = "This patient did not give his/her permission for reporting of pharmacogenomics results.";
            return TableUtil.createNoConsentReportTable(title,
                    noConsent,
                    TableUtil.TABLE_BOTTOM_MARGIN,
                    ReportResources.CONTENT_WIDTH_WIDE);
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