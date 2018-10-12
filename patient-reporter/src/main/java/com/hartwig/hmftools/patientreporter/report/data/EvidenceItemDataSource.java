package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRange;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRangeEvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public abstract class EvidenceItemDataSource {
    private static final Logger LOGGER = LogManager.getLogger(EvidenceItemDataSource.class);

    public static final FieldBuilder<?> GENE = field("Gene", String.class);
    public static final FieldBuilder<?> VARIANT = field("variant", String.class);
    public static final FieldBuilder<?> IMPACT = field("impact", String.class);
    public static final FieldBuilder<?> DRUG = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE = field("drugs type", String.class);
    public static final FieldBuilder<?> LEVEL = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE = field("response", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> LABEL = field("label", String.class);

    private EvidenceItemDataSource() {
    }

    @NotNull
    public static JRDataSource fromActionabilityVariants(@NotNull Map<EnrichedSomaticVariant, VariantEvidenceItems> evidenceItems,
            @NotNull Map<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> evidenceItemsRange) {
        final DRDataSource actionabilityVariantsDatasource = new DRDataSource(GENE.getName(),
                VARIANT.getName(),
                IMPACT.getName(),
                DRUG.getName(),
                DRUGS_TYPE.getName(),
                LEVEL.getName(),
                RESPONSE.getName(),
                SOURCE.getName(),
                LABEL.getName());

        for (Map.Entry<EnrichedSomaticVariant, VariantEvidenceItems> entry : evidenceItems.entrySet()) {

            String codingEffect = entry.getKey().canonicalHgvsCodingImpact();
            String proteinImpact = entry.getKey().canonicalHgvsProteinImpact();

            for (EvidenceItem variant : entry.getValue().onLabel()) {
                actionabilityVariantsDatasource.add(variant.gene(),
                        codingEffect,
                        proteinImpact,
                        variant.drug(),
                        variant.drugsType(),
                        variant.level(),
                        variant.response(),
                        sourceName(variant.source()),
                        "yes");
            }

            for (EvidenceItem variant : entry.getValue().offLabel()) {
                actionabilityVariantsDatasource.add(variant.gene(),
                        codingEffect,
                        proteinImpact,
                        variant.drug(),
                        variant.drugsType(),
                        variant.level(),
                        variant.response(),
                        sourceName(variant.source()),
                        "no");
            }
        }

        for (Map.Entry<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> entry : evidenceItemsRange.entrySet()) {

            String codingEffect = entry.getKey().canonicalHgvsCodingImpact();
            String proteinImpact = entry.getKey().canonicalHgvsProteinImpact();

            for (ActionabilityRange variantRange : entry.getValue().onLabel()) {
                actionabilityVariantsDatasource.add(variantRange.gene(),
                        codingEffect,
                        proteinImpact,
                        variantRange.drug(),
                        variantRange.drugsType(),
                        variantRange.level(),
                        variantRange.response(),
                        sourceName(variantRange.source()),
                        "yes");
            }

            for (ActionabilityRange variantRange : entry.getValue().offLabel()) {
                actionabilityVariantsDatasource.add(variantRange.gene(),
                        codingEffect,
                        proteinImpact,
                        variantRange.drug(),
                        variantRange.drugsType(),
                        variantRange.level(),
                        variantRange.response(),
                        sourceName(variantRange.source()),
                        "no");
            }
        }
        return actionabilityVariantsDatasource;
    }

    @NotNull
    private static String sourceName(@NotNull String source) {
        String sourceName = "";
        if (source.equals("oncoKb")) {
            sourceName = "OncoKB";
        } else if (source.equals("iclusion")) {
            sourceName = "Iclusion";
        } else if (source.equals("civic")) {
            sourceName = "CiViC";
        } else if (source.equals("cgi")) {
            sourceName = "CGI";
        }
        return sourceName;
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink(@NotNull List<EvidenceItem> evidenceItems) {
        List<String> linkReference = Lists.newArrayList();
        List<String> linkGene = Lists.newArrayList();
        for (EvidenceItem variant : evidenceItems) {
            linkReference.add(variant.reference());
            linkGene.add(variant.gene());
        }
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                final String source = data.getValue(SOURCE.getName());
                switch (source) {
                    case "oncoKb":
                        return "http://oncokb.org/#/gene/";
                                //+ linkGene.get(0) + "/alteration/" + linkReference.get(0);
                    case "iclusion":
                        return "https://www.iclusion.org";
                    case "cgi":
                        return "https://www.cancergenomeinterpreter.org/biomarkers";
                    case "civic":
                      //  String[] link = linkReference.get(0).split(":");
                        return "https://civic.genome.wustl.edu/links/variants/";
                    default:
                        return "";
                }
            }
        };
    }

    @NotNull
    public static FieldBuilder<?>[] actionabilityFields() {
        return new FieldBuilder<?>[] { GENE, VARIANT, IMPACT, DRUG, DRUGS_TYPE, LEVEL, RESPONSE, SOURCE, LABEL };
    }
}
