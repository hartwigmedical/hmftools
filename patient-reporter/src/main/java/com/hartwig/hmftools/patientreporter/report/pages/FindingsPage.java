package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.monospaceBaseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.report.components.COSMICLinkExpression;
import com.hartwig.hmftools.patientreporter.report.components.DataExpression;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.components.MicrosatelliteSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalLoadSection;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.VariantDataSource;

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
    abstract SequencedPatientReport report();

    @NotNull
    abstract HmfReporterData reporterData();

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(MainPageTopSection.buildWithImpliedPurity("HMF Sequencing Report",
                report().sampleReport(),
                report().impliedPurityString()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                pointMutationReport(report(), reporterData().drupFilter()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneCopyNumberReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneFusionReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                microsatelliteReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                mutationalLoadReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                geneDisruptionReport(report()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> pointMutationReport(@NotNull final SequencedPatientReport report,
            @NotNull final DrupFilter drupFilter) {
        final String geneMutationAddition = "Marked genes (*) are included in the DRUP study and indicate potential "
                + "eligibility in DRUP. Please note that the marking is NOT based on the specific mutation reported for "
                + "this sample, but only on a gene-level.";

        final ComponentBuilder<?, ?> table =
                !report.variants().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(VariantDataSource.variantFields())
                        .columns(col.column("Gene", VariantDataSource.GENE_FIELD),
                                col.column("Position", VariantDataSource.POSITION_FIELD),
                                col.column("Variant", VariantDataSource.VARIANT_FIELD),
                                col.column("Depth (VAF)", VariantDataSource.DEPTH_VAF_FIELD),
                                col.componentColumn("Predicted Effect", predictedEffectColumn()),
                                col.column("Cosmic", VariantDataSource.COSMIC_FIELD)
                                        .setHyperLink(hyperLink(new COSMICLinkExpression()))
                                        .setStyle(linkStyle()),
                                col.column("Ploidy (TAF)", VariantDataSource.PLOIDY_TAF_FIELD)))
                        .setDataSource(VariantDataSource.fromVariants(report.variants(), drupFilter))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table,
                cmp.verticalGap(15),
                cmp.horizontalList(cmp.horizontalGap(10),
                        cmp.text("*").setStyle(fontStyle()).setWidth(2),
                        cmp.text(geneMutationAddition).setStyle(fontStyle().setFontSize(8))));
    }

    @NotNull
    private static ComponentBuilder<?, ?> predictedEffectColumn() {
        return cmp.verticalList(cmp.horizontalList(cmp.text(DataExpression.fromField(VariantDataSource.HGVS_CODING_FIELD)),
                cmp.text(DataExpression.fromField(VariantDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(VariantDataSource.CONSEQUENCE_FIELD)));
    }

    @NotNull
    private static ComponentBuilder<?, ?> mutationalLoadReport(@NotNull SequencedPatientReport report) {
        return MutationalLoadSection.build(report.mutationalLoad());
    }

    @NotNull
    private static ComponentBuilder<?, ?> microsatelliteReport(@NotNull SequencedPatientReport report) {
        return MicrosatelliteSection.build(report.microsatelliteIndicator());
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneCopyNumberReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                !report.geneCopyNumbers().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(GeneCopyNumberDataSource.copyNumberFields())
                        .columns(col.column("Position", GeneCopyNumberDataSource.POSITION_FIELD),
                                col.column("Gene", GeneCopyNumberDataSource.GENE_FIELD),
                                col.column("Type", GeneCopyNumberDataSource.GAIN_OR_LOSS_FIELD),
                                col.column("Copies", GeneCopyNumberDataSource.COPY_NUMBER_FIELD))
                        .setDataSource(GeneCopyNumberDataSource.fromCopyNumbers(report.geneCopyNumbers())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gains & Losses").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneDisruptionReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.geneDisruptions().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(GeneDisruptionDataSource.geneDisruptionFields())
                .columns(col.column("Position", GeneDisruptionDataSource.POSITION_FIELD),
                        col.column("Gene", GeneDisruptionDataSource.GENE_FIELD),
                        col.column("Type", GeneDisruptionDataSource.TYPE_FIELD),
                        col.column("Context", GeneDisruptionDataSource.GENE_CONTEXT),
                        col.column("Copies", GeneDisruptionDataSource.COPIES))
                .setDataSource(GeneDisruptionDataSource.fromGeneDisruptions(report.geneDisruptions())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Disruptions").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> geneFusionReport(@NotNull final SequencedPatientReport report) {
        final ComponentBuilder<?, ?> table =
                !report.geneFusions().isEmpty()
                        ? cmp.subreport(monospaceBaseTable().fields(GeneFusionDataSource.geneFusionFields())
                        .columns(col.column("5' Gene", GeneFusionDataSource.GENE_FIELD),
                                col.column("5' Gene Context", GeneFusionDataSource.GENE_CONTEXT),
                                col.column("3' Gene", GeneFusionDataSource.PARTNER_GENE_FIELD),
                                col.column("3' Gene Context", GeneFusionDataSource.PARTNER_CONTEXT_FIELD),
                                col.column("Copies", GeneFusionDataSource.COPIES))
                        .setDataSource(GeneFusionDataSource.fromGeneFusions(report.geneFusions())))
                        : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Somatic Gene Fusions").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_DISTANCE),
                table);
    }
}
