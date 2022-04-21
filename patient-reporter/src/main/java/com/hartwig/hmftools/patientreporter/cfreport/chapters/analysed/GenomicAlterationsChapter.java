package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.GainsAndLosses;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneDisruptions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.HomozygousDisruptions;
import com.hartwig.hmftools.patientreporter.cfreport.data.Pharmacogenetics;
import com.hartwig.hmftools.patientreporter.cfreport.data.SomaticVariants;
import com.hartwig.hmftools.patientreporter.cfreport.data.TumorPurity;
import com.hartwig.hmftools.patientreporter.cfreport.data.ViralPresence;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.TextAlignment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class GenomicAlterationsChapter implements ReportChapter {

    // TODO Remove this toggle-off once we can remove position (blocked by DEV-810)
    private static final boolean DISPLAY_CLONAL_COLUMN = false;

    @NotNull
    private final AnalysedPatientReport patientReport;

    @NotNull
    private final SampleReport sampleReport;

    public GenomicAlterationsChapter(@NotNull final AnalysedPatientReport patientReport, @NotNull final SampleReport sampleReport) {
        this.patientReport = patientReport;
        this.sampleReport = sampleReport;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    public String name() {
        return "Genomic alteration details";
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        GenomicAnalysis genomicAnalysis = patientReport.genomicAnalysis();
        boolean hasReliablePurity = genomicAnalysis.hasReliablePurity();

        reportDocument.add(createPloidyPloidyTable(genomicAnalysis.averageTumorPloidy(),
                genomicAnalysis.impliedPurity(),
                hasReliablePurity));

        reportDocument.add(createTumorVariantsTable(genomicAnalysis.reportableVariants(),
                genomicAnalysis.notifyGermlineStatusPerVariant(),
                hasReliablePurity));

        reportDocument.add(createGainsAndLossesTable(genomicAnalysis.gainsAndLosses(),
                hasReliablePurity,
                genomicAnalysis.cnPerChromosome()));
        reportDocument.add(createFusionsTable(genomicAnalysis.geneFusions(), hasReliablePurity));
        reportDocument.add(createHomozygousDisruptionsTable(genomicAnalysis.homozygousDisruptions()));
        reportDocument.add(createDisruptionsTable(genomicAnalysis.geneDisruptions(), hasReliablePurity));
        reportDocument.add(createVirusTable(genomicAnalysis.reportableViruses(), sampleReport.reportViralPresence()));
        reportDocument.add(createPeachGenotypesTable(patientReport.peachGenotypes(), sampleReport.reportPharmogenetics()));
    }

    @NotNull
    private static Table createPloidyPloidyTable(double ploidy, double purity, boolean hasReliablePurity) {
        String title = "Tumor purity & ploidy";

        Table contentTable = TableUtil.createReportContentSmallTable(new float[] { 90, 95, 50 }, new Cell[] {});

        double impliedPurityPercentage = MathUtil.mapPercentage(purity, TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
        renderTumorPurity(hasReliablePurity,
                DataUtil.formatPercentage(impliedPurityPercentage),
                purity,
                TumorPurity.RANGE_MIN,
                TumorPurity.RANGE_MAX,
                contentTable);

        String copyNumber = GeneUtil.copyNumberToString(ploidy, hasReliablePurity);
        contentTable.addCell(TableUtil.createContentCell("Average tumor ploidy"));
        if (copyNumber.equals(DataUtil.NA_STRING)) {
            contentTable.addCell(TableUtil.createContentCell(copyNumber).setTextAlignment(TextAlignment.CENTER));
        } else {
            contentTable.addCell(TableUtil.createContentCellPurityPloidy(copyNumber).setTextAlignment(TextAlignment.CENTER));

        }
        contentTable.addCell(TableUtil.createContentCell(Strings.EMPTY));

        return TableUtil.createWrappingReportTable(title, null, contentTable);
    }

    @NotNull
    private static InlineBarChart createInlineBarChart(double value, double min, double max) {
        InlineBarChart chart = new InlineBarChart(value, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }

    private static void renderTumorPurity(boolean hasReliablePurity, @NotNull String valueLabel, double value, double min, double max,
            @NotNull Table table) {

        String label = "Tumor purity";
        table.addCell(TableUtil.createContentCell(label));

        if (hasReliablePurity) {
            table.addCell(TableUtil.createContentCellPurityPloidy(valueLabel).setTextAlignment(TextAlignment.CENTER));
            table.addCell(TableUtil.createContentCell(createInlineBarChart(value, min, max))
                    .setPadding(8)
                    .setTextAlignment(TextAlignment.CENTER));
        } else {
            table.addCell(TableUtil.createContentCell(Lims.PURITY_NOT_RELIABLE_STRING));
            table.addCell(TableUtil.createContentCell(Strings.EMPTY));
        }
    }

    @NotNull
    private static Table createTumorVariantsTable(@NotNull List<ReportableVariant> reportableVariants,
            @NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant, boolean hasReliablePurity) {
        String title = "Tumor specific variants";
        if (reportableVariants.isEmpty()) {
            return TableUtil.createNoneReportTable(title, null);
        }

        Table contentTable;
        if (DISPLAY_CLONAL_COLUMN) {
            contentTable = TableUtil.createReportContentTable(new float[] { 60, 70, 80, 70, 60, 40, 30, 60, 60, 50, 50 },
                    new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Position"),
                            TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("Protein"),
                            TableUtil.createHeaderCell("Read depth").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("tVAF").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Biallelic").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Hotspot").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Clonal").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Driver").setTextAlignment(TextAlignment.CENTER) });
        } else {
            contentTable = TableUtil.createReportContentTable(new float[] { 60, 70, 80, 70, 60, 40, 30, 60, 60, 50 },
                    new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Position"),
                            TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("Protein"),
                            TableUtil.createHeaderCell("Read depth").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("tVAF").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Biallelic").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Hotspot").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Driver").setTextAlignment(TextAlignment.CENTER) });
        }

        for (ReportableVariant variant : SomaticVariants.sort(reportableVariants)) {
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.geneDisplayString(variant,
                    notifyGermlineStatusPerVariant.get(variant))));
            contentTable.addCell(TableUtil.createContentCell(variant.gDNA()));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.determineCanonicalImpact(variant.canonicalHgvsCodingImpact(),
                    variant.canonicalEffect())));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.proteinAnnotationDisplayString(variant.canonicalHgvsProteinImpact(),
                    variant.canonicalEffect())));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(
                    variant.alleleReadCount() + " / ").setFont(ReportResources.fontBold())
                    .add(new Text(String.valueOf(variant.totalReadCount())).setFont(ReportResources.fontRegular()))
                    .setTextAlignment(TextAlignment.CENTER)));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.copyNumberString(variant.totalCopyNumber(), hasReliablePurity))
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.tVAFString(variant.tVAF(),
                    hasReliablePurity,
                    variant.totalCopyNumber())).setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.biallelicString(variant.biallelic(), hasReliablePurity))
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(SomaticVariants.hotspotString(variant.hotspot()))
                    .setTextAlignment(TextAlignment.CENTER));
            if (DISPLAY_CLONAL_COLUMN) {
                contentTable.addCell(TableUtil.createContentCell(SomaticVariants.clonalString(variant.clonalLikelihood()))
                        .setTextAlignment(TextAlignment.CENTER));
            }
            contentTable.addCell(TableUtil.createContentCell(variant.driverLikelihoodInterpretation().display()))
                    .setTextAlignment(TextAlignment.CENTER);
        }

        if (SomaticVariants.hasNotifiableGermlineVariant(notifyGermlineStatusPerVariant)) {
            contentTable.addCell(TableUtil.createLayoutCell(1, contentTable.getNumberOfColumns())
                    .add(new Paragraph("\n# Marked variant(s) are also present in the germline of the patient. Referral to a genetic "
                            + "specialist should be advised.").addStyle(ReportResources.subTextStyle())));
        }

        if (SomaticVariants.hasPhasedVariant(reportableVariants)) {
            contentTable.addCell(TableUtil.createLayoutCell(1, contentTable.getNumberOfColumns())
                    .add(new Paragraph("\n+ Marked protein (p.) annotation is based on multiple phased variants").addStyle(ReportResources.subTextStyle())));
        }

        return TableUtil.createWrappingReportTable(title, null, contentTable);
    }

    @NotNull
    private static Table createGainsAndLossesTable(@NotNull List<ReportableGainLoss> gainsAndLosses, boolean hasReliablePurity,
            @NotNull List<CnPerChromosomeArmData> cnPerChromosome) {
        String title = "Tumor specific gains & losses";
        if (gainsAndLosses.isEmpty()) {
            return TableUtil.createNoneReportTable(title, null);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 80, 80, 100, 80, 60, 60, 150 },
                new Cell[] { TableUtil.createHeaderCell("Chromosome"), TableUtil.createHeaderCell("Region"),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Type"), TableUtil.createHeaderCell("min copies"),
                        TableUtil.createHeaderCell("max copies"),
                        TableUtil.createHeaderCell("Chromosome arm copies").setTextAlignment(TextAlignment.CENTER) });

        List<ReportableGainLoss> sortedGainsAndLosses = GainsAndLosses.sort(gainsAndLosses);
        for (ReportableGainLoss gainLoss : sortedGainsAndLosses) {
            contentTable.addCell(TableUtil.createContentCell(gainLoss.chromosome()));
            contentTable.addCell(TableUtil.createContentCell(gainLoss.chromosomeBand()));
            contentTable.addCell(TableUtil.createContentCell(gainLoss.gene()));
            contentTable.addCell(TableUtil.createContentCell(gainLoss.interpretation().display()));
            contentTable.addCell(TableUtil.createContentCell(hasReliablePurity ? String.valueOf(gainLoss.minCopies()) : DataUtil.NA_STRING)
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(hasReliablePurity ? String.valueOf(gainLoss.maxCopies()) : DataUtil.NA_STRING)
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(GainsAndLosses.chromosomeArmCopyNumber(cnPerChromosome, gainLoss))
                    .setTextAlignment(TextAlignment.CENTER));
        }

        return TableUtil.createWrappingReportTable(title, null, contentTable);
    }

    @NotNull
    private static Table createHomozygousDisruptionsTable(@NotNull List<ReportableHomozygousDisruption> homozygousDisruptions) {
        String title = "Tumor specific homozygous disruptions";
        String subtitle = "Complete loss of wild type allele";
        if (homozygousDisruptions.isEmpty()) {
            return TableUtil.createNoneReportTable(title, subtitle);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 80, 80, 100 },
                new Cell[] { TableUtil.createHeaderCell("Chromosome"), TableUtil.createHeaderCell("Region"),
                        TableUtil.createHeaderCell("Gene"), });

        for (ReportableHomozygousDisruption homozygousDisruption : HomozygousDisruptions.sort(homozygousDisruptions)) {
            contentTable.addCell(TableUtil.createContentCell(homozygousDisruption.chromosome()));
            contentTable.addCell(TableUtil.createContentCell(homozygousDisruption.chromosomeBand()));
            contentTable.addCell(TableUtil.createContentCell(homozygousDisruption.gene()));
        }

        return TableUtil.createWrappingReportTable(title, subtitle, contentTable);
    }

    @NotNull
    private static Table createFusionsTable(@NotNull List<LinxFusion> fusions, boolean hasReliablePurity) {
        String title = "Tumor specific gene fusions";
        if (fusions.isEmpty()) {
            return TableUtil.createNoneReportTable(title, null);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 80, 80, 80, 40, 40, 40, 65, 70, 40 },
                new Cell[] { TableUtil.createHeaderCell("Fusion"), TableUtil.createHeaderCell("5' Transcript"),
                        TableUtil.createHeaderCell("3' Transcript"), TableUtil.createHeaderCell("5' End"),
                        TableUtil.createHeaderCell("3' Start"), TableUtil.createHeaderCell("Copies").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Phasing").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Type").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Driver").setTextAlignment(TextAlignment.CENTER), });

        for (LinxFusion fusion : GeneFusions.sort(fusions)) {
            contentTable.addCell(TableUtil.createContentCell(GeneFusions.name(fusion)));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(fusion.geneTranscriptStart()))
                    .addStyle(ReportResources.dataHighlightLinksStyle())
                    .setAction(PdfAction.createURI(GeneFusions.transcriptUrl(fusion.geneTranscriptStart()))));
            contentTable.addCell(TableUtil.createContentCell(new Paragraph(fusion.geneTranscriptEnd()).addStyle(ReportResources.dataHighlightLinksStyle())
                    .setAction(PdfAction.createURI(GeneFusions.transcriptUrl(fusion.geneTranscriptEnd())))));
            contentTable.addCell(TableUtil.createContentCell(fusion.geneContextStart()));
            contentTable.addCell(TableUtil.createContentCell(fusion.geneContextEnd()));
            contentTable.addCell(TableUtil.createContentCell(GeneUtil.copyNumberToString(fusion.junctionCopyNumber(), hasReliablePurity))
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(fusion.phased().displayStr()).setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(GeneFusions.displayStr(fusion.reportedType()))
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(fusion.likelihood().displayStr()).setTextAlignment(TextAlignment.CENTER));
        }

        return TableUtil.createWrappingReportTable(title, null, contentTable);
    }

    @NotNull
    private static Table createDisruptionsTable(@NotNull List<ReportableGeneDisruption> disruptions, boolean hasReliablePurity) {
        String title = "Tumor specific gene disruptions";
        if (disruptions.isEmpty()) {
            return TableUtil.createNoneReportTable(title, null);
        }

        Table contentTable = TableUtil.createReportContentTable(new float[] { 60, 50, 100, 50, 80, 85, 85 },
                new Cell[] { TableUtil.createHeaderCell("Location"), TableUtil.createHeaderCell("Gene"),
                        TableUtil.createHeaderCell("Disrupted range"),
                        TableUtil.createHeaderCell("Type").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Cluster ID").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Disrupted copies").setTextAlignment(TextAlignment.CENTER),
                        TableUtil.createHeaderCell("Undisrupted copies").setTextAlignment(TextAlignment.CENTER) });

        for (ReportableGeneDisruption disruption : GeneDisruptions.sort(disruptions)) {
            contentTable.addCell(TableUtil.createContentCell(disruption.location()));
            contentTable.addCell(TableUtil.createContentCell(disruption.gene()));
            contentTable.addCell(TableUtil.createContentCell(disruption.range()));
            contentTable.addCell(TableUtil.createContentCell(disruption.type())).setTextAlignment(TextAlignment.CENTER);
            contentTable.addCell(TableUtil.createContentCell(String.valueOf(disruption.clusterId()))
                    .setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(GeneUtil.copyNumberToString(disruption.junctionCopyNumber(),
                    hasReliablePurity)).setTextAlignment(TextAlignment.CENTER));
            contentTable.addCell(TableUtil.createContentCell(GeneUtil.copyNumberToString(disruption.undisruptedCopyNumber(),
                    hasReliablePurity)).setTextAlignment(TextAlignment.CENTER));
        }
        return TableUtil.createWrappingReportTable(title, null, contentTable);
    }

    @NotNull
    private static Table createVirusTable(@NotNull List<AnnotatedVirus> viruses, boolean reportViralPresence) {
        String title = "Tumor specific viral insertions";

        if (!reportViralPresence) {
            String noConsent = "This patient did not give his/her permission for reporting of virus results.";
            return TableUtil.createNoConsentReportTable(title, noConsent);
        } else if (viruses.isEmpty()) {
            return TableUtil.createNoneReportTable(title, null);
        } else {
            Table contentTable = TableUtil.createReportContentTable(new float[] { 150, 160, 100, 40 },
                    new Cell[] { TableUtil.createHeaderCell("Virus"),
                            TableUtil.createHeaderCell("Number of detected integration sites").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Viral coverage").setTextAlignment(TextAlignment.CENTER),
                            TableUtil.createHeaderCell("Driver").setTextAlignment(TextAlignment.CENTER) });

            for (AnnotatedVirus virus : viruses) {
                contentTable.addCell(TableUtil.createContentCell(virus.name()));
                contentTable.addCell(TableUtil.createContentCell(ViralPresence.createIntegrationSiteString(virus.integrations()))
                        .setTextAlignment(TextAlignment.CENTER));
                contentTable.addCell(TableUtil.createContentCell(ViralPresence.createViralCoverageString(virus.percentageCovered()))
                        .setTextAlignment(TextAlignment.CENTER));
                contentTable.addCell(TableUtil.createContentCell(virus.virusDriverLikelihoodType().display())
                        .setTextAlignment(TextAlignment.CENTER));
            }

            return TableUtil.createWrappingReportTable(title, null, contentTable);
        }
    }

    @NotNull
    private static Table createPeachGenotypesTable(@NotNull List<PeachGenotype> peachGenotypes, boolean reportPeach) {
        String title = "Pharmacogenetics";

        if (reportPeach) {
            if (peachGenotypes.isEmpty()) {
                return TableUtil.createNoneReportTable(title, null);
            } else {
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
                return TableUtil.createWrappingReportTable(title, null, contentTable);
            }
        } else {
            String noConsent = "This patient did not give his/her permission for reporting of pharmacogenomics results.";
            return TableUtil.createNoConsentReportTable(title, noConsent);
        }
    }
}