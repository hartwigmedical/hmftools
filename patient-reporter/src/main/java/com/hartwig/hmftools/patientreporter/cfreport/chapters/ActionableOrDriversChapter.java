package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.*;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
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

public class ActionableOrDriversChapter implements ReportChapter {

    private final AnalysedPatientReport patientReport;

    public ActionableOrDriversChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    public String getName() {
        return "Actionable or drivers";
    }

    public void render(@NotNull Document reportDocument) {

        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        reportDocument.add(createTumorVariantsTable("Tumor specific variants",
                patientReport.somaticVariants(), hasReliablePurityFit));
        reportDocument.add(createGainsAndLossesTable("Tumor specific gains & losses",
                patientReport.geneCopyNumbers(), hasReliablePurityFit));
        reportDocument.add(createSomaticFusionsTable( "Somatic gene fusions",
                patientReport.geneFusions(), hasReliablePurityFit));
        reportDocument.add(createDisruptionsTable("Tumor specific gene disruptions",
                patientReport.geneDisruptions(), hasReliablePurityFit));
    }

    @NotNull
    private Table createTumorVariantsTable(@NotNull final String title, @NotNull final List<ReportableSomaticVariant> somaticVariants, boolean hasReliablePurityFit) {

        final List<ReportableSomaticVariant> sortedVariants = SomaticVariants.sort(somaticVariants);
        if (sortedVariants.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
        }

        // Create content table
        Table contentTable = TableUtil.createReportContentTable(new float[] {45, 75, 50, 60, 40, 60, 40, 50, 50, 35}, new Cell[]  {
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

            InlineBarChart chart = new InlineBarChart(hasReliablePurityFit ? variant.adjustedVAF() : 0, 0, 1);
            chart.setEnabled(hasReliablePurityFit);
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
            contentTable.addCell(TableUtil.getContentCell("Medium")); //SomaticVariants.getDriverString(variant.driverLikelihood())
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
            contentTable.addCell(TableUtil.getContentCell(hasReliablePurityFit ? String.valueOf(copyNumber.value()) : DataUtil.NAString).setTextAlignment(TextAlignment.RIGHT));
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
            contentTable.addCell(TableUtil.getContentCell(GeneUtil.getPloidyToCopiesString(fusion.ploidy(), hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell("")); // Spacer
            contentTable.addCell(TableUtil.getContentCell(new Paragraph(fusion.source())
                    .setAction(PdfAction.createURI(GeneFusions.getSourceUrl(fusion.source())))));
        }


        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

    @NotNull
    private Table createDisruptionsTable(@NotNull final String title, @NotNull final List<ReportableGeneDisruption> disruptions, boolean hasReliablePurityFit) {

        final List<ReportableGeneDisruption> sortedDisruptions = GeneDisruptions.sort(disruptions);
        if (sortedDisruptions.isEmpty()) {
            return TableUtil.createNoneReportTable(title);
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

        for (ReportableGeneDisruption disruption : sortedDisruptions) {
            contentTable.addCell(TableUtil.getContentCell(disruption.location()));
            contentTable.addCell(TableUtil.getContentCell(disruption.gene()));
            contentTable.addCell(TableUtil.getContentCell(disruption.range()));
            contentTable.addCell(TableUtil.getContentCell(disruption.type()));
            contentTable.addCell(TableUtil.getContentCell(GeneUtil.getPloidyToCopiesString(disruption.ploidy(), hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell(GeneDisruptions.getCopyNumberString(disruption.geneMinCopies(), hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableUtil.getContentCell(GeneDisruptions.getCopyNumberString(disruption.geneMaxCopies(), hasReliablePurityFit)).setTextAlignment(TextAlignment.RIGHT));
        }

        // Create report table that handles page breaks
        return TableUtil.createWrappingReportTable(title, contentTable);

    }

}
