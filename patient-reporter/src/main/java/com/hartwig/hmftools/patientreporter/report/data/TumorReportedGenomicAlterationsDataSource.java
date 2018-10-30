package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class TumorReportedGenomicAlterationsDataSource {
    public static final FieldBuilder<?> NUMBER_SOMATIC_VARIANTS = field("Number of somatic variants", String.class);
    public static final FieldBuilder<?> GENES_SOMATIC_VARIANTS_WITH_DRIVER_LIKELIHOOD =
            field("Somatic variants with driver likelihood > 0.8", String.class);
    public static final FieldBuilder<?> GENES_WITH_AMPLIFICATION = field("Genes with amplification", String.class);
    public static final FieldBuilder<?> GENES_WITH_LOSS = field("Genes with (partial) loss", String.class);
    public static final FieldBuilder<?> FUSION_GENES = field("Fusion genes", String.class);
    public static final FieldBuilder<?> DISRUPTION_GENES = field("Disruption genes", String.class);


    private TumorReportedGenomicAlterationsDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] TumorReportedGenomicAlterationsSummaryFields() {
        return new FieldBuilder<?>[] { NUMBER_SOMATIC_VARIANTS, GENES_SOMATIC_VARIANTS_WITH_DRIVER_LIKELIHOOD, GENES_WITH_AMPLIFICATION,
                GENES_WITH_LOSS, FUSION_GENES, DISRUPTION_GENES };
    }

    @NotNull
    public static JRDataSource fromTumorReportedGenomicAlterationsSummary(@NotNull AnalysedPatientReport report) {
        final DRDataSource TumorReportedGenomicAlterationsSummaryDataSource = new DRDataSource(NUMBER_SOMATIC_VARIANTS.getName(),
                GENES_SOMATIC_VARIANTS_WITH_DRIVER_LIKELIHOOD.getName(),
                GENES_WITH_AMPLIFICATION.getName(),
                GENES_WITH_LOSS.getName(),
                FUSION_GENES.getName(),
                DISRUPTION_GENES.getName());

        TumorReportedGenomicAlterationsSummaryDataSource.add(countSomaticVariants(report.somaticVariants()),
                somaticVariantsWithDriver(report.somaticVariants()),
                amplificationGenes(report.geneCopyNumbers()),
                lossGenes(report.geneCopyNumbers()),
                geneFusions(report.geneFusions()),
                geneDisruptions(report.geneDisruptions()));

        return TumorReportedGenomicAlterationsSummaryDataSource;
    }

    @NotNull
    private static String countSomaticVariants(@NotNull List<ReportableSomaticVariant> variants) {
        int countSomaticVariants = 0;
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            countSomaticVariants++;
        }
        return Integer.toString(countSomaticVariants);
    }

    @NotNull
    private static String somaticVariantsWithDriver(@NotNull List<ReportableSomaticVariant> variants) {
        List<String> somaticVariants = Lists.newArrayList();
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            if (SomaticVariantDataSource.driverField(variant).equals("High")) {
                somaticVariants.add(variant.gene());
            }
        }
        return String.join("\n", somaticVariants);
    }

    @NotNull
    private static String amplificationGenes(@NotNull final List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersAmplification = Lists.newArrayList();
        for (final GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("gain")) {
                geneCopyNumbersAmplification.add(copyNumber.gene());
            }
        }
        return String.join("\n", geneCopyNumbersAmplification);
    }

    @NotNull
    private static String lossGenes(@NotNull final List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersLoss = Lists.newArrayList();
        for (final GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("full loss") || GeneCopyNumberDataSource.type(copyNumber)
                    .equals("partial loss")) {
                geneCopyNumbersLoss.add(copyNumber.gene());
            }
        }
        return String.join("\n", geneCopyNumbersLoss);
    }

    @NotNull
    private static String geneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        List<String> geneFusions = Lists.newArrayList();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            geneFusions.add(GeneFusionDataSource.name(fusion));
        }
        return String.join("\n", geneFusions);
    }

    @NotNull
    private static String geneDisruptions(@NotNull List<ReportableGeneDisruption> disruptions) {
        List<String> geneDisruptions = Lists.newArrayList();
        for (ReportableGeneDisruption disruption : GeneDisruptionDataSource.sort(disruptions)) {
            geneDisruptions.add(disruption.gene());
        }
        return String.join("\n", geneDisruptions);
    }
}

