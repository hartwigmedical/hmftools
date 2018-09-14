package com.hartwig.hmftools.patientreporter.report.data;

import static com.google.common.base.Strings.repeat;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class SomaticVariantDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_DETAILS_FIELD = field("variant_details", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> IS_HOTSPOT_FIELD = field("is_hotspot", String.class);
    public static final FieldBuilder<?> PLOIDY_VAF_FIELD = field("ploidy_vaf", String.class);
    public static final FieldBuilder<?> CLONAL_PERCENTAGE_FIELD = field("clonal_probability", String.class);
    public static final FieldBuilder<?> WILDTYPE_STATUS_FIELD = field("wildtype_status", String.class);
    public static final FieldBuilder<?> DRIVER_PROBABILITY_FIELD = field("driver_probability", String.class);
    public static final FieldBuilder<?> ACTIONABILITY_LEVEL_FIELD = field("actionability_level", String.class);

    private SomaticVariantDataSource() {
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull final List<EnrichedSomaticVariant> variants,
            @NotNull DrupFilter drupFilter) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_DETAILS_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                IS_HOTSPOT_FIELD.getName(),
                PLOIDY_VAF_FIELD.getName(),
                CLONAL_PERCENTAGE_FIELD.getName(),
                WILDTYPE_STATUS_FIELD.getName(),
                DRIVER_PROBABILITY_FIELD.getName(),
                ACTIONABILITY_LEVEL_FIELD.getName());

        for (final EnrichedSomaticVariant variant : variants) {
            final String displayGene = drupFilter.test(variant) ? variant.gene() + " *" : variant.gene();

            variantDataSource.add(displayGene,
                    variantDetailsField(variant),
                    readDepthField(variant),
                    variant.hotspot() ? "Yes" : "No",
                    PatientReportFormat.correctValueForFitStatus(fitStatus, ploidyVafField(variant)),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, PatientReportFormat.formatPercentWithDefaultCutoffs(0D)),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, Strings.EMPTY),
                    PatientReportFormat.formatPercentWithDefaultCutoffs(0D),
                    Strings.EMPTY);
        }

        return variantDataSource;
    }

    @NotNull
    private static String variantDetailsField(@NotNull EnrichedSomaticVariant variant) {
        return variant.canonicalHgvsCodingImpact() + " (" + variant.canonicalHgvsProteinImpact() + ")";
    }

    @NotNull
    private static String readDepthField(@NotNull EnrichedSomaticVariant variant) {
        return variant.alleleReadCount() + " / " + variant.totalReadCount() + " ("
                + PatientReportFormat.formatPercent(variant.alleleFrequency()) + ")";
    }

    @NotNull
    private static String ploidyVafField(@NotNull EnrichedSomaticVariant variant) {
        return descriptiveBAF(variant.adjustedCopyNumber(), variant.minorAllelePloidy()) + " ("
                + PatientReportFormat.formatPercent(variant.adjustedVAF()) + ")";
    }

    @NotNull
    @VisibleForTesting
    static String descriptiveBAF(double adjustedCopyNumber, double minorAllelePloidy) {
        int totalAlleleCount = (int) Math.max(0, Math.round(adjustedCopyNumber));
        int minorAlleleCount = (int) Math.max(0, Math.round(minorAllelePloidy));
        int majorAlleleCount = Math.max(0, totalAlleleCount - minorAlleleCount);

        return formatBafField("A", Math.max(minorAlleleCount, majorAlleleCount)) + formatBafField("B",
                Math.min(minorAlleleCount, majorAlleleCount));
    }

    @NotNull
    private static String formatBafField(@NotNull final String allele, final int count) {
        return count < 10 ? repeat(allele, count) : allele + "[" + count + "x]";
    }

    @NotNull
    public static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_DETAILS_FIELD, READ_DEPTH_FIELD, IS_HOTSPOT_FIELD, PLOIDY_VAF_FIELD,
                CLONAL_PERCENTAGE_FIELD, WILDTYPE_STATUS_FIELD, DRIVER_PROBABILITY_FIELD, ACTIONABILITY_LEVEL_FIELD };
    }
}
