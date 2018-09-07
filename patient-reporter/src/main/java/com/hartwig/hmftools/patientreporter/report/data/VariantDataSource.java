package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public class VariantDataSource {

    private static final String COSMIC_IDENTIFIER = "COSM";

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> VARIANT_TYPE_FIELD = field("variant_type", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> KNOWLEDGEBASE_KEY_FIELD = field("knowledgebase_key", String.class);
    private static final FieldBuilder<?> KNOWLEDGEBASE_URL_FIELD = field("knowledgebase_url", String.class);
    public static final FieldBuilder<?> PLOIDY_FIELD = field("ploidy", String.class);
    public static final FieldBuilder<?> VAF_FIELD = field("vaf", String.class);
    public static final FieldBuilder<?> CLONAL_FIELD = field("clonal", String.class);
    public static final FieldBuilder<?> WILDTYPE_FIELD = field("wildtype", String.class);
    public static final FieldBuilder<?> DRIVER_FIELD = field("driver", String.class);
    public static final FieldBuilder<?> ACTIONABLE_FIELD = field("actionable", String.class);

    private VariantDataSource() {
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull final List<VariantReport> variantReports,
            @NotNull DrupFilter drupFilter) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_FIELD.getName(),
                VARIANT_TYPE_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                KNOWLEDGEBASE_KEY_FIELD.getName(),
                KNOWLEDGEBASE_URL_FIELD.getName(),
                PLOIDY_FIELD.getName(),
                VAF_FIELD.getName(),
                CLONAL_FIELD.getName(),
                WILDTYPE_FIELD.getName(),
                DRIVER_FIELD.getName(),
                ACTIONABLE_FIELD.getName());

        for (final VariantReport variantReport : variantReports) {
            final String displayGene = drupFilter.test(variantReport) ? variantReport.gene() + " *" : variantReport.gene();
            variantDataSource.add(displayGene,
                    variantReport.proteinImpact(),
                    variantReport.proteinImpactType(),
                    variantReport.readDepthField(),
                    variantReport.knowledgebaseKey(),
                    cosmicUrl(variantReport.knowledgebaseKey()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.ploidy()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.purityAdjustedVAFField()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.clonalityStatus()),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, variantReport.wildTypeStatus()),
                    variantReport.driverStatus(),
                    variantReport.actionabilityStatus());
        }

        return variantDataSource;
    }

    @NotNull
    private static String cosmicUrl(@Nullable final String cosmicID) {
        if (cosmicID == null) {
            return Strings.EMPTY;
        }
        final int identifierPos = cosmicID.indexOf(COSMIC_IDENTIFIER);
        String cosmicIdentifier;
        if (identifierPos >= 0) {
            cosmicIdentifier = cosmicID.substring(identifierPos + COSMIC_IDENTIFIER.length());
        } else {
            cosmicIdentifier = cosmicID;
        }

        return "http://cancer.sanger.ac.uk/cosmic/mutation/overview?genome=37&id=" + cosmicIdentifier;
    }

    @NotNull
    public static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_FIELD, VARIANT_TYPE_FIELD, READ_DEPTH_FIELD, KNOWLEDGEBASE_KEY_FIELD,
                KNOWLEDGEBASE_URL_FIELD, PLOIDY_FIELD, VAF_FIELD, CLONAL_FIELD, WILDTYPE_FIELD, DRIVER_FIELD, ACTIONABLE_FIELD };
    }

    @NotNull
    public static AbstractSimpleExpression<String> knowledgeBaseHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return data.getValue(KNOWLEDGEBASE_URL_FIELD.getName());
            }
        };
    }

}
