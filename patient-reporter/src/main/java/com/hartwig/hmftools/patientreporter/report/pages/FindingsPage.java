package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.monospaceBaseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import java.util.List;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.actionability.ReportableClinicalTrials;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItems;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.ChordSection;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.components.MicrosatelliteSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalLoadSection;
import com.hartwig.hmftools.patientreporter.report.components.TumorMutationBurdenSection;
import com.hartwig.hmftools.patientreporter.report.data.ClinicalTrialDataSource;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceItemDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GermlineVariantDataSource;
import com.hartwig.hmftools.patientreporter.report.data.SomaticVariantDataSource;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class FindingsPage {

    private static final int HEADER_TO_TABLE_DISTANCE = 6;

    @NotNull
    abstract AnalysedPatientReport report();

    @NotNull
    abstract SequencedReportData reporterData();

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(MainPageTopSection.buildWithImpliedPurity(Commons.TITLE_SEQUENCE,
                report().sampleReport(),
                impliedPurityString(report())),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                evidenceItemReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                clinicalTrialToReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                somaticVariantReport(report(), reporterData().panelGeneModel()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                germlineVariantReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneCopyNumberReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneFusionReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                microsatelliteReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                Double.toString(report().chordValue().iterator().next().hrdValue()) != null ? chordReport(report()) : cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                tumorMutationalLoadReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                tumorMutationalBurdenReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneDisruptionReport(report()));
    }

    @NotNull
    private static String impliedPurityString(@NotNull AnalysedPatientReport report) {
        return report.fitStatus() == FittedPurityStatus.NO_TUMOR
                ? "[below detection threshold]"
                : PatientReportFormat.formatPercent(report.impliedPurity());
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceItemReport(@NotNull AnalysedPatientReport report) {
        List<EvidenceItem> reportableItems = ReportableEvidenceItems.reportableEvidenceItems(report.evidenceItems());

        final ComponentBuilder<?, ?> table =
                reportableItems.size() > 0
                        ? cmp.subreport(monospaceBaseTable().fields(EvidenceItemDataSource.evidenceItemFields())
                        .columns(col.column("Event", EvidenceItemDataSource.EVENT_FIELD),
                                col.column("Drug", EvidenceItemDataSource.DRUG_FIELD),
                                col.column("Drugs type", EvidenceItemDataSource.DRUGS_TYPE_FIELD),
                                col.column("Level", EvidenceItemDataSource.LEVEL_FIELD),
                                col.column("Response", EvidenceItemDataSource.RESPONSE_FIELD),
                                col.column("Source", EvidenceItemDataSource.SOURCE_FIELD)
                                        .setHyperLink(hyperLink(EvidenceItemDataSource.sourceHyperlink()))
                                        .setStyle(linkStyle()),
                                col.column("On-Label", EvidenceItemDataSource.ON_LABEL_FIELD))
                        .setDataSource(EvidenceItemDataSource.fromEvidenceItems(reportableItems)))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Clinical Evidence").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> clinicalTrialToReport(@NotNull AnalysedPatientReport report) {
        List<EvidenceItem> clinicalTrials = ReportableClinicalTrials.reportableTrials(report.evidenceItems());

        final ComponentBuilder<?, ?> table =
                clinicalTrials.size() > 0
                        ? cmp.subreport(monospaceBaseTable().fields(ClinicalTrialDataSource.clinicalTrialFields())
                        .columns(col.column("Event", ClinicalTrialDataSource.EVENT_FIELD),
                                col.column("Trial", ClinicalTrialDataSource.TRIAL_FIELD),
                                col.column("Source", ClinicalTrialDataSource.SOURCE_FIELD)
                                        .setHyperLink(hyperLink(ClinicalTrialDataSource.sourceHyperlink()))
                                        .setStyle(linkStyle()),
                                col.column("CCMO", ClinicalTrialDataSource.CCMO_FIELD),
                                col.column("On-Label", ClinicalTrialDataSource.ON_LABEL_FIELD))
                        .setDataSource(ClinicalTrialDataSource.fromClinicalTrials(clinicalTrials)))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Clinical Trials").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> somaticVariantReport(@NotNull AnalysedPatientReport report, @NotNull GeneModel panelGeneModel) {
        final String drupEligibilityAddition = "Marked genes (*) are included in the DRUP study and indicate potential "
                + "eligibility in DRUP. Please note that the marking is NOT based on the specific mutation reported for "
                + "this sample, but only on a gene-level.";

        final ComponentBuilder<?, ?> table =
                !report.somaticVariants().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(SomaticVariantDataSource.fields())
                        .columns(col.column("Gene", SomaticVariantDataSource.GENE_FIELD),
                                col.column("Variant", SomaticVariantDataSource.VARIANT_FIELD).setFixedWidth(90),
                                col.column("Impact", SomaticVariantDataSource.IMPACT_FIELD).setFixedWidth(80),
                                col.column("Read Depth", SomaticVariantDataSource.READ_DEPTH_FIELD),
                                col.column("Hotspot", SomaticVariantDataSource.IS_HOTSPOT_FIELD),
                                col.column("Ploidy (VAF)", SomaticVariantDataSource.PLOIDY_VAF_FIELD).setFixedWidth(80),
                                col.column("Clonality", SomaticVariantDataSource.CLONAL_STATUS_FIELD),
                                col.column("Biallelic", SomaticVariantDataSource.BIALLELIC_FIELD),
                                col.column("Driver", SomaticVariantDataSource.DRIVER_FIELD)))
                        .setDataSource(SomaticVariantDataSource.fromVariants(report.fitStatus(),
                                report.somaticVariants(),
                                report.somaticVariantDriverCatalog(),
                                panelGeneModel))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table,
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("*").setStyle(fontStyle()).setWidth(2),
                        cmp.text(drupEligibilityAddition).setStyle(fontStyle().setFontSize(8))));
    }

    @NotNull
    private static ComponentBuilder<?, ?> germlineVariantReport(@NotNull AnalysedPatientReport report) {
        String noVariantsFoundText = report.hasGermlineAnalysis() ? "None" : "Germline analysis is not available";
        final ComponentBuilder<?, ?> table =
                !report.germlineVariants().isEmpty() && report.hasGermlineAnalysis()
                        ? cmp.subreport(monospaceBaseTable().fields(GermlineVariantDataSource.fields())
                        .columns(col.column("Gene", GermlineVariantDataSource.GENE_FIELD),
                                col.column("Variant", GermlineVariantDataSource.VARIANT_FIELD).setFixedWidth(90),
                                col.column("Impact", GermlineVariantDataSource.IMPACT_FIELD).setFixedWidth(80),
                                col.column("Read Depth", GermlineVariantDataSource.READ_DEPTH_FIELD),
                                col.column("Germline Status", GermlineVariantDataSource.GERMLINE_STATUS_FIELD),
                                col.column("Ploidy (VAF)", GermlineVariantDataSource.PLOIDY_VAF_FIELD).setFixedWidth(80),
                                col.column("Biallelic", GermlineVariantDataSource.BIALLELIC_FIELD))
                        .setDataSource(GermlineVariantDataSource.fromVariants(report.fitStatus(), report.germlineVariants())))
                        : cmp.text(noVariantsFoundText).setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Actionable Germline Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneCopyNumberReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                !report.geneCopyNumbers().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(GeneCopyNumberDataSource.copyNumberFields())
                        .columns(col.column("Chromosome", GeneCopyNumberDataSource.CHROMOSOME),
                                col.column("Chromosome band", GeneCopyNumberDataSource.CHROMOSOME_BAND),
                                col.column("Gene", GeneCopyNumberDataSource.GENE_FIELD),
                                col.column("Type", GeneCopyNumberDataSource.GAIN_OR_LOSS_FIELD),
                                col.column("Copies", GeneCopyNumberDataSource.COPY_NUMBER_FIELD))
                        .setDataSource(GeneCopyNumberDataSource.fromCopyNumbers(report.fitStatus(), report.geneCopyNumbers())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gains & Losses").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                !report.geneFusions().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(GeneFusionDataSource.geneFusionFields())
                        .columns(col.column("Fusion", GeneFusionDataSource.FUSION_FIELD),
                                col.column("5' Transcript", GeneFusionDataSource.START_TRANSCRIPT_FIELD)
                                        .setHyperLink(hyperLink(GeneFusionDataSource.transcriptUrl(GeneFusionDataSource.START_TRANSCRIPT_FIELD)))
                                        .setStyle(linkStyle()),
                                col.column("3' Transcript", GeneFusionDataSource.END_TRANSCRIPT_FIELD)
                                        .setHyperLink(hyperLink(GeneFusionDataSource.transcriptUrl(GeneFusionDataSource.END_TRANSCRIPT_FIELD)))
                                        .setStyle(linkStyle()),
                                col.column("5' End", GeneFusionDataSource.START_CONTEXT_FIELD),
                                col.column("3' Start", GeneFusionDataSource.END_CONTEXT_FIELD),
                                col.column("Copies", GeneFusionDataSource.COPIES_FIELD),
                                col.column("Source", GeneFusionDataSource.SOURCE_FIELD)
                                        .setHyperLink(hyperLink(GeneFusionDataSource.sourceHyperlink()))
                                        .setStyle(linkStyle()))
                        .setDataSource(GeneFusionDataSource.fromGeneFusions(report.fitStatus(), report.geneFusions())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Fusions").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> microsatelliteReport(@NotNull AnalysedPatientReport report) {
        return MicrosatelliteSection.build(report.microsatelliteIndelsPerMb(), report.fitStatus());
    }

    @NotNull
    private static ComponentBuilder<?, ?> chordReport(@NotNull AnalysedPatientReport report) {
        return ChordSection.build(report.chordValue().iterator().next().hrdValue());
    }

    @NotNull
    private static ComponentBuilder<?, ?> tumorMutationalLoadReport(@NotNull AnalysedPatientReport report) {
        return MutationalLoadSection.build(report.tumorMutationalLoad(), report.fitStatus());
    }

    @NotNull
    private static ComponentBuilder<?, ?> tumorMutationalBurdenReport(@NotNull AnalysedPatientReport report) {
        return TumorMutationBurdenSection.build(report.tumorMutationalBurden(), report.fitStatus());
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.geneDisruptions().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(GeneDisruptionDataSource.geneDisruptionFields())
                .columns(col.column("Location", GeneDisruptionDataSource.LOCATION_FIELD),
                        col.column("Gene", GeneDisruptionDataSource.GENE_FIELD),
                        col.column("Disrupted Range", GeneDisruptionDataSource.RANGE_FIELD).setFixedWidth(120),
                        col.column("Type", GeneDisruptionDataSource.TYPE_FIELD),
                        col.column("Copies", GeneDisruptionDataSource.COPIES_FIELD),
                        col.column("Gene Min Copies", GeneDisruptionDataSource.GENE_MIN_COPIES).setFixedWidth(80),
                        col.column("Gene Max Copies", GeneDisruptionDataSource.GENE_MAX_COPIES).setFixedWidth(80))
                .setDataSource(GeneDisruptionDataSource.fromGeneDisruptions(report.fitStatus(), report.geneDisruptions())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Disruptions").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }
}