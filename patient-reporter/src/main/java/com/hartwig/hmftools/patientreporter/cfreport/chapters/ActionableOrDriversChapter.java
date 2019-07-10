package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneCopyNumbers;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneDisruptions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.SomaticVariants;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

public class ActionableOrDriversChapter implements ReportChapter {

    @NotNull
    private final AnalysedPatientReport patientReport;

    public ActionableOrDriversChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @Override
    @NotNull
    public String name() {
        return "Genomic alteration details";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        reportDocument.add(createTumorVariantsTable(patientReport.reportableVariants(), hasReliablePurityFit));
        reportDocument.add(createGainsAndLossesTable(patientReport.geneCopyNumbers(), hasReliablePurityFit));
        reportDocument.add(createSomaticFusionsTable(patientReport.geneFusions(), hasReliablePurityFit));
        reportDocument.add(createDisruptionsTable(patientReport.geneDisruptions(), hasReliablePurityFit));
    }

    @NotNull
    private static Table createTumorVariantsTable(@NotNull List<ReportableVariant> reportableVariants, boolean hasReliablePurityFit) {
        final String title = "Tumor specific variants";
        if (reportableVariants.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 50, 75, 60, 60, 45, 65, 35, 50, 45, 35 },
                new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("gDNA"), TableUtil.createHeaderCell("Variant"),
                        TableUtil.createHeaderCell("Impact"),
                        TableUtil.createHeaderCell("Read depth").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Hotspot"), TableUtil.createHeaderCell("Ploidy (VAF)"), TableUtil.createHeaderCell(),
                        TableUtil.createHeaderCell("Clonality"), TableUtil.createHeaderCell("Biallelic"),
                        TableUtil.createHeaderCell("Driver") });

        final List<ReportableVariant> sortedVariants = SomaticVariants.sort(reportableVariants);
        for (ReportableVariant variant : sortedVariants) {
            InlineBarChart chart = new InlineBarChart(hasReliablePurityFit ? variant.adjustedVAF() : 0, 0, 1);
            chart.enabled(hasReliablePurityFit);
            chart.setWidth(20);
            chart.setHeight(4);

            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.geneDisplayString(variant)));
            contentTable.addCell(TableUtil.createContentCell(variant.gDNA()));
            contentTable.addCell(TableUtil.createContentCell(variant.hgvsCodingImpact()));
            contentTable.addCell(TableUtil.createContentCell(variant.hgvsProteinImpact()));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(
                    variant.alleleReadCount() + " / ").setFont(ReportResources.fontBold())
                    .add(new Text(String.valueOf(variant.totalReadCount())).setFont(ReportResources.fontRegular()))
                    .setTextAlignment(TextAlignment.CENTER)));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.hotspotString(variant.hotspot())));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.ploidyVaf(variant.adjustedCopyNumber(),
                    variant.minorAllelePloidy(),
                    variant.adjustedVAF(),
                    hasReliablePurityFit)));
            contentTable.addCell(TableUtil.createContentCell(chart).setVerticalAlignment(VerticalAlignment.MIDDLE));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.clonalityString(variant.clonality(), hasReliablePurityFit)));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.biallelicString(variant.biallelic(),
                    variant.driverCategory(),
                    hasReliablePurityFit)));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.driverString(variant.driverLikelihood())));
        }

        contentTable.addCell(TableUtil.createLayoutCell(1, contentTable.getNumberOfColumns())
                .setPaddingTop(10)
                .add(new Paragraph("* Marked gene(s) are included in the DRUP study and indicate potential eligibility in "
                        + "DRUP. Please note that the marking is NOT based on the specific mutation").addStyle(ReportResources.subTextStyle()))
                .add(new Paragraph("reported for this sample, but only on a gene-level.").setPaddingLeft(5)
                        .addStyle(ReportResources.subTextStyle())));

        if (SomaticVariants.hasNotifiableGermlineVariant(reportableVariants)) {
            contentTable.addCell(TableUtil.createLayoutCell(1, contentTable.getNumberOfColumns())
                    .add(new Paragraph("# Marked variant(s) are also present in the germline of the patient. "
                            + "Referral to a genetic specialist should be advised.").addStyle(ReportResources.subTextStyle())));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }

    @NotNull
    private static Table createGainsAndLossesTable(@NotNull List<GeneCopyNumber> copyNumbers, boolean hasReliablePurityFit) {
        final String title = "Tumor specific gains & losses";
        if (copyNumbers.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 80, 100, 80, 45, 125 },
                new Cell[] { TableUtil.createHeaderCell("Chromosome"), TableUtil.createHeaderCell("Chromosome band"),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Type"),
                        TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT), TableUtil.createHeaderCell("") });

        final List<GeneCopyNumber> sortedCopyNumbers = GeneCopyNumbers.sort(copyNumbers);
        for (GeneCopyNumber copyNumber : sortedCopyNumbers) {
            Long reportableCopyNumber = Math.round(Math.max(0, copyNumber.minCopyNumber()));
            contentTable.addCell(TableUtil.createContentCell(copyNumber.chromosome()));
            contentTable.addCell(TableUtil.createContentCell(copyNumber.chromosomeBand()));
            contentTable.addCell(TableUtil.createContentCell(copyNumber.gene()));
            contentTable.addCell(TableUtil.createContentCell(GeneCopyNumbers.type(copyNumber)));
            contentTable.addCell(TableUtil.createContentCell(hasReliablePurityFit
                    ? String.valueOf(reportableCopyNumber)
                    : DataUtil.NA_STRING).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.createContentCell(""));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }

    @NotNull
    private static Table createSomaticFusionsTable(@NotNull List<ReportableGeneFusion> fusions, boolean hasReliablePurityFit) {
        final String title = "Somatic gene fusions";
        if (fusions.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 90, 82.5f, 82.5f, 37.5f, 37.5f, 40, 30, 100 },
                new Cell[] { TableUtil.createHeaderCell("Fusion"), TableUtil.createHeaderCell("5' Transcript"),
                        TableUtil.createHeaderCell("3' Transcript"), TableUtil.createHeaderCell("5' End"),
                        TableUtil.createHeaderCell("3' Start"), TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                        TableUtil.createHeaderCell(""), TableUtil.createHeaderCell("Source") });

        final List<ReportableGeneFusion> sortedFusions = GeneFusions.sort(fusions);
        for (ReportableGeneFusion fusion : sortedFusions) {
            contentTable.addCell(TableUtil.createContentCell(GeneFusions.name(fusion)));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(fusion.geneTranscriptStart()))
                    .addStyle(ReportResources.dataHighlightLinksStyle())
                    .setAction(PdfAction.createURI(GeneFusions.transcriptUrl(fusion.geneTranscriptStart()))));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(fusion.geneTranscriptEnd()).addStyle(ReportResources.dataHighlightLinksStyle())
                    .setAction(PdfAction.createURI(GeneFusions.transcriptUrl(fusion.geneTranscriptEnd())))));
            contentTable.addCell(TableUtil.createContentCell(fusion.geneContextStart()));
            contentTable.addCell(TableUtil.createContentCell(fusion.geneContextEnd()));
            contentTable.addCell(TableUtil.createContentCell(GeneUtil.ploidyToCopiesString(fusion.ploidy(), hasReliablePurityFit))
                    .setTextAlignment(TextAlignment.RIGHT));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }

    @NotNull
    private static Table createDisruptionsTable(@NotNull List<ReportableGeneDisruption> disruptions, boolean hasReliablePurityFit) {
        final String title = "Tumor specific gene disruptions";
        if (disruptions.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 80, 100, 80, 40, 65, 65 },
                new Cell[] { TableUtil.createHeaderCell("Location"), TableUtil.createHeaderCell("Gene"),
                        TableUtil.createHeaderCell("Disrupted range"), TableUtil.createHeaderCell("Type"),
                        TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                        TableUtil.createHeaderCell("Gene \nmin copies").setTextAlignment(TextAlignment.RIGHT),
                        TableUtil.createHeaderCell("Gene \nmax copies").setTextAlignment(TextAlignment.RIGHT) });

        final List<ReportableGeneDisruption> sortedDisruptions = GeneDisruptions.sort(disruptions);
        for (ReportableGeneDisruption disruption : sortedDisruptions) {
            contentTable.addCell(TableUtil.createContentCell(disruption.location()));
            contentTable.addCell(TableUtil.createContentCell(disruption.gene()));
            contentTable.addCell(TableUtil.createContentCell(disruption.range()));
            contentTable.addCell(TableUtil.createContentCell(disruption.type()));
            contentTable.addCell(TableUtil.createContentCell(GeneUtil.ploidyToCopiesString(disruption.ploidy(), hasReliablePurityFit))
                    .setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.createContentCell(GeneDisruptions.copyNumberString(disruption.geneMinCopies(),
                    hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.createContentCell(GeneDisruptions.copyNumberString(disruption.geneMaxCopies(),
                    hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
        }

        return TableUtil.createWrappingReportTable(title, contentTable);
    }
}
