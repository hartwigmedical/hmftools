package com.hartwig.hmftools.patientreporter.report.data;

import static com.google.common.base.Strings.repeat;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.algo.DriverProbabilityModel;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class SomaticVariantDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> IMPACT_FIELD = field("impact", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> IS_HOTSPOT_FIELD = field("is_hotspot", String.class);
    public static final FieldBuilder<?> PLOIDY_VAF_FIELD = field("ploidy_vaf", String.class);
    public static final FieldBuilder<?> CLONAL_STATUS_FIELD = field("clonal_status", String.class);
    public static final FieldBuilder<?> BIALLELIC_FIELD = field("biallelic", String.class);
    public static final FieldBuilder<?> DRIVER_PROBABILITY_FIELD = field("driver_probability", String.class);

    private static final double MIN_PERCENTAGE_CUTOFF_DRIVER_PROB = 0.2;
    private static final double MAX_PERCENTAGE_CUTOFF_DRIVER_PROB = 0.9;

    private SomaticVariantDataSource() {
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull List<EnrichedSomaticVariant> variants,
            @NotNull DriverProbabilityModel driverProbabilityModel, @NotNull GeneModel panelGeneModel) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_FIELD.getName(),
                IMPACT_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                IS_HOTSPOT_FIELD.getName(),
                PLOIDY_VAF_FIELD.getName(),
                CLONAL_STATUS_FIELD.getName(),
                BIALLELIC_FIELD.getName(),
                DRIVER_PROBABILITY_FIELD.getName());

        for (final EnrichedSomaticVariant variant : variants) {
            final DriverCategory driverCategory = panelGeneModel.geneDriverCategory(variant.gene());

            final String displayGene =
                    panelGeneModel.drupActionableGenes().contains(variant.gene()) ? variant.gene() + " *" : variant.gene();

            String biallelic = Strings.EMPTY;
            if (driverCategory != null && driverCategory == DriverCategory.TSG) {
                biallelic = variant.biallelic() ? "Yes" : "No";
            }

            DriverCatalog driver = driverProbabilityModel.catalogForVariant(variant);
            String driverProbabilityString = driver != null ? PatientReportFormat.formatPercentWithCutoffs(driver.driverLikelihood(),
                    MIN_PERCENTAGE_CUTOFF_DRIVER_PROB,
                    MAX_PERCENTAGE_CUTOFF_DRIVER_PROB) : "N/A";

            variantDataSource.add(displayGene,
                    variant.canonicalHgvsCodingImpact(),
                    variant.canonicalHgvsProteinImpact(),
                    readDepthField(variant),
                    driverCategory != null && driverCategory == DriverCategory.ONCO ? hotspotField(variant) : Strings.EMPTY,
                    PatientReportFormat.correctValueForFitStatus(fitStatus, ploidyVafField(variant)),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, clonalityField(variant)),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, biallelic),
                    driverProbabilityString);
        }

        return variantDataSource;
    }

    @NotNull
    private static String hotspotField(@NotNull SomaticVariant variant) {
        switch (variant.hotspot()) {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            case NON_HOTSPOT:
                return "No";
            default:
                return Strings.EMPTY;
        }
    }

    @NotNull
    private static String clonalityField(@NotNull EnrichedSomaticVariant variant) {
        switch (variant.clonality()) {
            case CLONAL:
                return "Clonal";
            case SUBCLONAL:
                return "Subclonal";
            case INCONSISTENT:
                return "Inconsistent";
            default:
                return Strings.EMPTY;
        }
    }

    @NotNull
    private static String readDepthField(@NotNull EnrichedSomaticVariant variant) {
        return variant.alleleReadCount() + " / " + variant.totalReadCount();
    }

    @NotNull
    private static String ploidyVafField(@NotNull EnrichedSomaticVariant variant) {
        return descriptiveBAF(variant.adjustedCopyNumber(), variant.minorAllelePloidy()) + " ("
                + PatientReportFormat.formatPercent(Math.min(1, variant.adjustedVAF())) + ")";
    }

    @NotNull
    @VisibleForTesting
    static String descriptiveBAF(double adjustedCopyNumber, double minorAllelePloidy) {
        int totalAlleleCount = (int) Math.max(0, Math.round(adjustedCopyNumber));
        int minorAlleleCount = (int) Math.max(0, Math.round(minorAllelePloidy));
        int majorAlleleCount = Math.max(0, totalAlleleCount - minorAlleleCount);

        return formatBAFField("A", Math.max(minorAlleleCount, majorAlleleCount)) + formatBAFField("B",
                Math.min(minorAlleleCount, majorAlleleCount));
    }

    @NotNull
    private static String formatBAFField(@NotNull String allele, int count) {
        return count < 10 ? repeat(allele, count) : allele + "[" + count + "x]";
    }

    @NotNull
    public static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_FIELD, READ_DEPTH_FIELD, IS_HOTSPOT_FIELD, PLOIDY_VAF_FIELD, CLONAL_STATUS_FIELD,
                BIALLELIC_FIELD, DRIVER_PROBABILITY_FIELD };
    }
}
