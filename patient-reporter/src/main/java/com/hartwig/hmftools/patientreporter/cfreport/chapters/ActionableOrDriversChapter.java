package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.*;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class ActionableOrDriversChapter extends ReportChapter {

    @Override
    public String getName() {
        return "Actionable or drivers";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(@NotNull final AnalysedPatientReport patientReport, @NotNull Document reportDocument) {

        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        reportDocument.add(createTumorVariantsTable("Tumor specific variants",
                patientReport.somaticVariants(), hasReliablePurityFit));
        reportDocument.add(createGainsAndLossesTable("Tumor specific gains & losses",
                patientReport.geneCopyNumbers(), hasReliablePurityFit));
        reportDocument.add(createSomaticFusionsTable( "Somatic gene fusions",
                patientReport.geneFusions(), hasReliablePurityFit));
        reportDocument.add(createDisruptionsTable());
    }

    @NotNull
    private Table createTumorVariantsTable(@NotNull final String title, @NotNull final List<ReportableSomaticVariant> somaticVariants, boolean hasReliablePurityFit) {

        final List<ReportableSomaticVariant> sortedVariants = SomaticVariants.sort(somaticVariants);
        if (sortedVariants.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {45, 75, 50, 60, 40, 60, 40, 55, 50, 30}, new Cell[]  {
                TableUtil.getHeaderCell("Gene"),
                TableUtil.getHeaderCell("Variant"),
                TableUtil.getHeaderCell("Impact"),
                TableUtil.getHeaderCell("Read depth").setTextAlignment(TextAlignment.CENTER),
                TableUtil.getHeaderCell("Hotspot"),
                TableUtil.getHeaderCell("Ploidy (VAF)"),
                TableUtil.getHeaderCell(), // Spacer for graph
                TableUtil.getHeaderCell("Clonality"),
                TableUtil.getHeaderCell("Biallelic"),
                TableUtil.getHeaderCell("Driver")
        });

        for (ReportableSomaticVariant variant : sortedVariants) {

            // @TODO: handle purity fit miss in chart
            InlineBarChart chart = new InlineBarChart(hasReliablePurityFit ? variant.adjustedVAF() : 0, 0, 1);
            chart.setWidth(20);
            chart.setHeight(4);

            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getGeneDisplayString(variant)));
            contentTable.addCell(TableUtil.getContentCell(variant.hgvsCodingImpact()));
            contentTable.addCell(TableUtil.getContentCell(variant.hgvsProteinImpact()));
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(variant.alleleReadCount() + " / ")
                    .setFont(ReportResources.getFontBold())
                    .add(new Text(String.valueOf(variant.totalReadCount()))
                            .setFont(ReportResources.getFontRegular()))
                    .setTextAlignment(TextAlignment.CENTER)));
            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getHotspotString(variant.hotspot())));
            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getPloidyVaf(variant.adjustedCopyNumber(), variant.minorAllelePloidy(), variant.adjustedVAF(), hasReliablePurityFit)));
            contentTable.addCell(TableUtil.getContentCell(chart).setVerticalAlignment(VerticalAlignment.MIDDLE));
            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getClonalityString(variant.clonality(), hasReliablePurityFit)));
            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getBiallelicString(variant.biallelic(), variant.driverCategory(), hasReliablePurityFit)));
            contentTable.addCell(TableUtil.getContentCell(SomaticVariants.getDriverString(variant.driverLikelihood())));
        }

        // Add table footnotes
        contentTable.addCell(TableUtil.getLayoutCell(1, contentTable.getNumberOfColumns())
                .setPaddingTop(10)
                .add(new Paragraph("* Marked gene(s) are included in the DRUP study and indicate potential eligibility in " +
                        "DRUP. Please note that the marking is NOT based on the specific mutation reported for this sample, " +
                        "but only on a gene-level.")
                        .addStyle(ReportResources.subTextStyle())));
        contentTable.addCell(TableUtil.getLayoutCell(1, contentTable.getNumberOfColumns())
                .add(new Paragraph("# Marked variant(s) are also present in the germline of the patient. Referral " +
                        "to a genetic specialist should be considered if a hereditary condition is suspected.")
                        .addStyle(ReportResources.subTextStyle())));

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private Table createGainsAndLossesTable(@NotNull final String title, @NotNull final List<GeneCopyNumber> copyNumbers, boolean hasReliablePurityFit) {

        final List<GeneCopyNumber> sortedCopyNumbers = GeneCopyNumbers.sort(copyNumbers);
        if (sortedCopyNumbers.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {60, 80, 100, 80, 45, 125}, new Cell[]  {
                TableUtil.getHeaderCell("Chromosome"),
                TableUtil.getHeaderCell("Chromosome band"),
                TableUtil.getHeaderCell("Gene"),
                TableUtil.getHeaderCell("Type"),
                TableUtil.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                TableUtil.getHeaderCell("") // Spacer
        });

        for (GeneCopyNumber copyNumber : copyNumbers) {
            contentTable.addCell(TableUtil.getContentCell(copyNumber.chromosome()));
            contentTable.addCell(TableUtil.getContentCell(copyNumber.chromosomeBand()));
            contentTable.addCell(TableUtil.getContentCell(copyNumber.gene()));
            contentTable.addCell(TableUtil.getContentCell(GeneCopyNumbers.getType(copyNumber)));
            contentTable.addCell(TableUtil.getContentCell(hasReliablePurityFit ? String.valueOf(copyNumber.value()) : Util.NAString).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell("")); // Spacer
        }

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private Table createSomaticFusionsTable(@NotNull final String title, @NotNull final List<ReportableGeneFusion> fusions, boolean hasReliablePurityFit) {

        final List<ReportableGeneFusion> sortedFusions = GeneFusions.sort(fusions);
        if (sortedFusions.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {90, 82.5f, 82.5f, 37.5f, 37.5f, 40, 30, 100}, new Cell[]  {
                TableUtil.getHeaderCell("Fusion"),
                TableUtil.getHeaderCell("5' Transcript"),
                TableUtil.getHeaderCell("3' Transcript"),
                TableUtil.getHeaderCell("5' End"),
                TableUtil.getHeaderCell("3' Start"),
                TableUtil.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                TableUtil.getHeaderCell(""), // Spacer
                TableUtil.getHeaderCell("Source")
        });

        for (ReportableGeneFusion fusion : sortedFusions) {
            contentTable.addCell(TableUtil.getContentCell(GeneFusions.getName(fusion)));
            contentTable.addCell(TableUtil.getContentCell(fusion.geneStartTranscript()));
            contentTable.addCell(TableUtil.getContentCell(fusion.geneEndTranscript()));
            contentTable.addCell(TableUtil.getContentCell(fusion.geneContextStart()));
            contentTable.addCell(TableUtil.getContentCell(fusion.geneContextEnd()));
            contentTable.addCell(TableUtil.getContentCell(GeneFusions.getPloidyToCopiesString(fusion.ploidy(), hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell("")); // Spacer
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(fusion.source())
                    .setAction(PdfAction.createURI(GeneFusions.getSourceUrl(fusion.source())))));
        }


        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private Table createDisruptionsTable() {

        final String chapterTitle = "Tumor specific gene disruptions";
        final boolean isAvailable = true;

        if (!isAvailable) {
            return TableUtil.createNoneReportTable(chapterTitle);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {60, 80, 100, 80, 40, 65, 65}, new Cell[]  {
                TableUtil.getHeaderCell("Location"),
                TableUtil.getHeaderCell("Gene"),
                TableUtil.getHeaderCell("Disrupted range"),
                TableUtil.getHeaderCell("Type"),
                TableUtil.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                TableUtil.getHeaderCell("Gene \nmin copies").setTextAlignment(TextAlignment.RIGHT),
                TableUtil.getHeaderCell("Gene \nmax copies").setTextAlignment(TextAlignment.RIGHT)
        });

        for (int i = 0; i < 4; i++) {
            contentTable.addCell(TableUtil.getContentCell("q23.31"));
            contentTable.addCell(TableUtil.getContentCell("PTEN"));
            contentTable.addCell(TableUtil.getContentCell("Intron 5 -> Intron 6"));
            contentTable.addCell(TableUtil.getContentCell("DEL"));
            contentTable.addCell(TableUtil.getContentCell("1.8").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell("2").setTextAlignment(TextAlignment.RIGHT));
        }

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(chapterTitle, contentTable);

    }

}
