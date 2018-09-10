package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public class VariantDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_DETAILS_FIELD = field("variant_details", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> COSMIC_ID_FIELD = field("cosmic_id", String.class);
    private static final FieldBuilder<?> COSMIC_URL_FIELD = field("cosmic_url", String.class);
    public static final FieldBuilder<?> PLOIDY_VAF_FIELD = field("ploidy_vaf", String.class);
    public static final FieldBuilder<?> CLONAL_PERCENTAGE_FIELD = field("clonal_probability", String.class);
    public static final FieldBuilder<?> WILDTYPE_STATUS_FIELD = field("wildtype_status", String.class);
    public static final FieldBuilder<?> DRIVER_PROBABILITY_FIELD = field("driver_probability", String.class);
    public static final FieldBuilder<?> ACTIONABILITY_LEVEL_FIELD = field("actionability_level", String.class);

    private VariantDataSource() {
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull final List<VariantReport> variantReports,
            @NotNull DrupFilter drupFilter) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_DETAILS_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                COSMIC_ID_FIELD.getName(),
                COSMIC_URL_FIELD.getName(),
                PLOIDY_VAF_FIELD.getName(),
                CLONAL_PERCENTAGE_FIELD.getName(),
                WILDTYPE_STATUS_FIELD.getName(),
                DRIVER_PROBABILITY_FIELD.getName(),
                ACTIONABILITY_LEVEL_FIELD.getName());

        for (final VariantReport variantReport : variantReports) {
            final String displayGene = drupFilter.test(variantReport) ? variantReport.gene() + " *" : variantReport.gene();
            variantDataSource.add(displayGene,
                    variantReport.variantDetails(),
                    variantReport.readDepth(),
                    variantReport.cosmicID(),
                    variantReport.cosmicUrl(),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.ploidyVaf()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus,
                            PatientReportFormat.formatPercent(variantReport.clonalProbability())),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.wildTypeStatus()),
                    PatientReportFormat.formatPercent(variantReport.driverProbability()),
                    variantReport.actionabilityLevel());
        }

        return variantDataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_DETAILS_FIELD, READ_DEPTH_FIELD, COSMIC_ID_FIELD, COSMIC_URL_FIELD, PLOIDY_VAF_FIELD,
                CLONAL_PERCENTAGE_FIELD, WILDTYPE_STATUS_FIELD, DRIVER_PROBABILITY_FIELD, ACTIONABILITY_LEVEL_FIELD };
    }

    @NotNull
    public static AbstractSimpleExpression<String> cosmicHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return data.getValue(COSMIC_URL_FIELD.getName());
            }
        };
    }
}
