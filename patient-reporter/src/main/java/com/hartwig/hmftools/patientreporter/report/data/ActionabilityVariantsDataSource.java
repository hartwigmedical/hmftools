package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public class ActionabilityVariantsDataSource {
    public static final FieldBuilder<?> EVENT = field("event", String.class);
    public static final FieldBuilder<?> CHROMOSOME = field("chromosome", String.class);
    public static final FieldBuilder<?> REF = field("ref", String.class);
    public static final FieldBuilder<?> ALT = field("alt", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> DRUG = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE = field("drugs type", String.class);
    public static final FieldBuilder<?> LEVEL = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE = field("response", String.class);
    public static final FieldBuilder<?> LABEL = field("label", String.class);

    private ActionabilityVariantsDataSource() {
    }

    @NotNull
    public static JRDataSource fromActionabilityVariants(@NotNull List<EvidenceItem> evidenceItems,
            @NotNull List<VariantEvidenceItems> label) {
        final DRDataSource actionabilityVariantsDatasource = new DRDataSource(EVENT.getName(),
                CHROMOSOME.getName(),
                REF.getName(),
                ALT.getName(),
                DRUG.getName(),
                DRUGS_TYPE.getName(),
                LEVEL.getName(),
                RESPONSE.getName(),
                LABEL.getName(),
                SOURCE.getName());

        for (EvidenceItem variant : evidenceItems) {
            actionabilityVariantsDatasource.add(variant.gene(),
                    variant.chromosome(),
                    variant.ref(),
                    variant.alt(),
                    variant.drug(),
                    variant.drugsType(),
                    variant.level(),
                    variant.response(),
                    "label",
                    variant.source());
        }
        return actionabilityVariantsDatasource;
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                final String source = data.getValue(SOURCE.getName());
                switch (source) {
                    case "oncoKb":
                        return "http://oncokb.org/#/";
                    case "iclusion":
                        return "https://http://www.iclusion.com";
                    case "cgi":
                        return "https://www.cancergenomeinterpreter.org/biomarkers";
                    case "civic":
                        return "https://civicdb.org/browse/somaticVariants";
                    default:
                        return "";
                }
            }
        };
    }

    @NotNull
    public static FieldBuilder<?>[] actionabilityFields() {
        return new FieldBuilder<?>[] { EVENT, CHROMOSOME, REF, ALT, DRUG, DRUGS_TYPE, LEVEL, RESPONSE, SOURCE, LABEL };
    }
}
