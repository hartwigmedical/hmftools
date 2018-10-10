package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityVariant;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class ActionabilityVariantsDataSource {
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> DRUG = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE = field("drugs type", String.class);
    public static final FieldBuilder<?> LEVEL = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE = field("response", String.class);

    private ActionabilityVariantsDataSource() {
    }

    @NotNull
    public static JRDataSource fromActionabilityVariants(@NotNull List<ActionabilityVariant> actionabilityVariants) {
        final DRDataSource actionabilityVariantsDatasource = new DRDataSource(
                SOURCE.getName(),
                DRUG.getName(),
                DRUGS_TYPE.getName(),
                LEVEL.getName(),
                RESPONSE.getName());

        for (ActionabilityVariant variant : actionabilityVariants) {
            actionabilityVariantsDatasource.add(variant.gene(), variant.source(), variant.drug(), variant.drugsType(), variant.level(), variant.response());
        }

        return actionabilityVariantsDatasource;
    }

    @NotNull
    public static FieldBuilder<?>[] actionabilityFields() {
        return new FieldBuilder<?>[] {SOURCE, DRUG, DRUGS_TYPE, LEVEL, RESPONSE };
    }
}
