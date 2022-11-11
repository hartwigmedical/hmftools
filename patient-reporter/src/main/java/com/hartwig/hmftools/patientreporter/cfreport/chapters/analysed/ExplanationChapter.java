package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ExplanationChapter implements ReportChapter {

    public ExplanationChapter() {
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @NotNull
    @Override
    public String name() {
        return "Report explanation";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 10, 1, 10, 1, 10, }));
        table.setWidth(contentWidth());

        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on the report in general")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on the reported clinical evidence")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported somatic variants")));

        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] { "The analysis is based on reference genome version GRCh37." }))
                .add(createContentDivWithLinkThree("The gene transcripts used for reporting can be downloaded from ",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip",
                        ". In general the used transcripts are the canonical transcripts as defined by Ensembl.",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip"))
                .add(createContentDiv(new String[] {
                        "Variant detection in samples with lower tumor content is less sensitive. In case of a low tumor "
                                + "purity (below 20%) likelihood of failing to detect potential variants increases.",
                        "The (implied) tumor purity is the percentage of tumor cells in the tumor material based on analysis of "
                                + "whole genome data." })));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDivWithLinkFour("The Clinical Knowledgebase (CKB) ",
                        "(https://ckbhome.jax.org/)",
                        " is used to annotate "
                                + "variants of all types with clinical evidence, with a hyperlink to the specific evidence items when "
                                + "available. The evidence is gathered from CKB without further checks or interpretation. "
                                + "This also means that if a certain evidence item or drug- biomarker is missing from the knowledgebase "
                                + "it will also not be included in this report.\n More details about CKB can be found in their Glossary "
                                + "Of Terms.",
                        "(https://ckbhome.jax.org/about/glossaryOfTerms)",
                        "" + "https://ckbhome.jax.org/",
                        "https://ckbhome.jax.org/about/glossaryOfTerms"))
                .add(createContentDivWithLinkThree("Clinical trials are matched against the iClusion database ",
                        "https://iclusion.org",
                        " including a link to the specific trial.\n",
                        "https://iclusion.org"))
                .add(createContentDiv(new String[] {
                        "Hartwig Medical Foundation is not responsible for the content of the " + "knowledgebases used to generate this "
                                + "report. Furthermore, Hartwig Medical Foundation is not liable and cannot be held accountable for any "
                                + "incorrectness, incompleteness or error of any other kind in the knowledgebases, or the external "
                                + "software used to harmonize and curate the knowledgebases." })));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] {
                        "The 'Read Depth' displays the raw number of reads supporting the variant versus the total "
                                + "number of reads on the mutated position.",
                        "The 'Copies' field indicates the number of alleles present in the tumor on this particular mutated position.",
                        "The 'tVAF' field displays the variant allele frequency corrected for tumor purity.",
                        "The 'Biallelic' field indicates whether the variant is present across all alleles in the tumor "
                                + "(and is including variants with loss-of-heterozygosity).",
                        "The 'Driver' field represents the driver probability on gene level and is calculated based on the HMF database. A "
                                + "variant in a gene with High driver likelihood is likely to be positively selected  "
                                + "during the oncogenic process." })));

        table.addCell(TableUtil.createLayoutCell(1, 5).setHeight(30));

        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported gene copy numbers")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported gene fusions")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported gene disruptions")));

        table.addCell(TableUtil.createLayoutCell()

                .add(createContentDiv(new String[] { "The lowest copy number value along the exonic regions of the canonical transcript is"
                        + " determined as a measure for the gene's copy number.",
                        "Copy numbers are corrected for the implied tumor purity and represent the number of copies in the tumor DNA.",
                        "Any gene with less than 0.5 copies along the entire canonical transcript is reported as a full loss.",
                        "Any gene where only a part along the canonical transcript has less than 0.5 copies is reported "
                                + "as a partial loss.",
                        "Any gene with more copies than 3 times the average tumor ploidy along the entire canonical transcript is reporte "
                                + " as a full gain.",
                        "Any gene where only a part of the canonical transcript has more copies than 3 times the average tumor ploidy "
                                + "is reported as a partial gain.", })));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] {
                        "The canonical, or otherwise longest transcript that is validly fused, is reported. " }))
                .add(createContentDivWithLinkThree("Fusions are restricted to a selection of known fusions and can be downloaded from ",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip ",
                        ".",
                        "https://storage.googleapis.com/hmf-public/OncoAct-Resources/latest_oncoact.zip"))
                .add(createContentDiv(new String[] {
                        "We additionally select fusions where one partner is promiscuous in either 5' or 3' position.",
                        "The 'Driver' field is set to HIGH in case the fusion is a known pathogenic fusion, or otherwise a fusion where "
                                + "the promiscuous partner is fused in an exon range that is typically observed in literature. \n"
                                + "All other fusions get assigned a LOW driver likelihood." })));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] {
                        "Genes are reported as being disrupted if their canonical transcript has been disrupted.",
                        "The range of the disruption is indicated by the intron/exon/promoter region of the break point "
                                + "and the direction the disruption faces.",
                        "The type of disruption can be INV (inversion), DEL (deletion), DUP (duplication), INS "
                                + "(insertion), SGL (single) or BND (translocation).",
                        "A gene for which no wild type exists anymore in the tumor DNA due to disruption(s) "
                                + "is reported in a separate section called 'homozygous disruptions'." })));

        // Is needed to set details on new page
        table.addCell(TableUtil.createLayoutCell(1, 5).setHeight(30));
        table.addCell(TableUtil.createLayoutCell(1, 5).setHeight(30));

        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported viral insertions")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported pharmacogenetics")));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createSectionTitle("Details on reported HLA Alleles")));

        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] { "Virusses will be reported if they are present in our reporting database as "
                        + "clinically relevant (HPV, MCV, HBV, EBV and HHV-8) and DNA integration for the virus can be detected. "
                        + "If the virus is clinically relevant and no DNA integration is found, the following conditions must be met:\n"
                        + "- Percentage covered of the viral genome is >90%\n"
                        + "- Coverage of the virus DNA is higher than expected tumor mean coverage\n",
                        "Reporting of EBV is independent of tumor integration. This means that to be reportable, the viral EBV"
                                + " genome must be covered >90% and the coverage of the virus must be higher than the expected clonal "
                                + "mean coverage." })));
        table.addCell(TableUtil.createLayoutCell());
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
        table.addCell(TableUtil.createLayoutCell()
                .add(createContentDiv(new String[] { "HLA Class I types (HLA-A, HLA-B and HLA-C) are "
                        + "reported based on blood analysis, and also includes the tumor status of each of those alleles (somatic mutations, "
                        + "complete loss, and/or allelic imbalance)\n" })

                        .add(createContentDivWithLinkThree("The IMGT/HLA database ",
                                "https://www.ebi.ac.uk/ipd/imgt/hla/",
                                " is used as a "
                                        + "reference set of Human MHC class I alleles. HLA typing is done to 4-digits, which means it uniquely "
                                        + "identifies a specific protein, but ignores synonymous variants (6 digits) and intronic differences "
                                        + "(8 digits).",
                                "https://www.ebi.ac.uk/ipd/imgt/hla/"))));

        reportDocument.add(table);
    }

    @NotNull
    private static Paragraph createSectionTitle(@NotNull String sectionTitle) {
        return new Paragraph(sectionTitle).addStyle(ReportResources.smallBodyHeadingStyle());
    }

    @NotNull
    private static Paragraph createParaGraphWithLinkFour(@NotNull String string1, @NotNull String string2, @NotNull String string3,
            @NotNull String string4, @NotNull String link1, @NotNull String link2) {
        return new Paragraph(string1).addStyle(ReportResources.subTextStyle())
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string2).addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(link1)))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string3).addStyle(ReportResources.subTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING)
                .add(new Text(string4).addStyle(ReportResources.urlStyle()).setAction(PdfAction.createURI(link2)))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
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

    @NotNull
    private static Div createContentDivWithLinkFour(@NotNull String string1, @NotNull String string2, @NotNull String string3,
            @NotNull String string4, @NotNull String link1, @NotNull String link2) {
        Div div = new Div();
        div.add(createParaGraphWithLinkFour(string1, string2, string3, string4, link1, link2));
        return div;
    }

    @NotNull
    private static Div createContentDiv(@NotNull String[] contentParagraphs) {
        Div div = new Div();
        for (String s : contentParagraphs) {
            div.add(new Paragraph(s).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING));
        }
        return div;
    }
}
