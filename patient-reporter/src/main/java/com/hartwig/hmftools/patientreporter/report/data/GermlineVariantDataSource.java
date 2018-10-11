package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GermlineVariantDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> IMPACT_FIELD = field("impact", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> GERMLINE_STATUS_FIELD = field("germline_status", String.class);
    public static final FieldBuilder<?> PLOIDY_VAF_FIELD = field("ploidy_vaf", String.class);
    public static final FieldBuilder<?> BIALLELIC_FIELD = field("biallelic", String.class);

    private GermlineVariantDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] fields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_FIELD, READ_DEPTH_FIELD, GERMLINE_STATUS_FIELD, PLOIDY_VAF_FIELD,
                BIALLELIC_FIELD };
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull List<GermlineVariant> variants) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_FIELD.getName(),
                IMPACT_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                GERMLINE_STATUS_FIELD.getName(),
                PLOIDY_VAF_FIELD.getName(),
                BIALLELIC_FIELD.getName());

        for (GermlineVariant variant : variants) {
            variantDataSource.add(variant.gene(),
                    variant.variant(),
                    variant.impact(),
                    variant.readDepth(),
                    variant.germlineStatus(),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variant.ploidyVaf()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variant.biallelic()));
        }

        return variantDataSource;
    }



}
