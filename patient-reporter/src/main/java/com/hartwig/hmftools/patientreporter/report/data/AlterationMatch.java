package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.apiclients.civic.data.CivicVariantType;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariantWithEvidence;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.definition.ReportParameters;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AlterationMatch {
    public static final FieldBuilder<?> MATCH_TYPE = field("matchType", String.class);
    public static final FieldBuilder<?> NAME = field("name", String.class);
    public static final FieldBuilder<?> CHROMOSOME = field("chromosome", String.class);
    public static final FieldBuilder<?> START = field("start", String.class);
    public static final FieldBuilder<?> STOP = field("stop", String.class);
    public static final FieldBuilder<?> CHROMOSOME2 = field("chromosome2", String.class);
    public static final FieldBuilder<?> START2 = field("start2", String.class);
    public static final FieldBuilder<?> STOP2 = field("stop2", String.class);
    public static final FieldBuilder<?> VARIANT_TYPE = field("variantType", String.class);
    public static final FieldBuilder<?> HGVS_EXPRESSIONS = field("hgvsExpressions", String.class);
    public static final FieldBuilder<?> SUMMARY_URL = field("summaryUrl", String.class);

    public abstract String getMatchType();

    @Nullable
    public abstract String getName();

    @Nullable
    public abstract String getChromosome();

    @Nullable
    public abstract String getStart();

    @Nullable
    public abstract String getStop();

    @Nullable
    public abstract String getChromosome2();

    @Nullable
    public abstract String getStart2();

    @Nullable
    public abstract String getStop2();

    public abstract String getVariantType();

    public abstract String getHgvsExpressions();

    public abstract String getSummaryUrl();

    public static AlterationMatch of(@NotNull final String matchType, @NotNull final CivicVariantWithEvidence variant) {
        return ImmutableAlterationMatch.of(matchType,
                variant.name(),
                variant.coordinates().chromosome(),
                Objects.toString(variant.coordinates().start(), ""),
                Objects.toString(variant.coordinates().stop(), ""),
                variant.coordinates().chromosome2(),
                Objects.toString(variant.coordinates().start2(), ""),
                Objects.toString(variant.coordinates().stop2(), ""),
                Strings.join(variant.variantTypes().stream().map(CivicVariantType::name).collect(Collectors.toList()), ','),
                Strings.join(variant.hgvsExpressions(), ','),
                variant.summaryUrl());
    }

    @NotNull
    public static AbstractSimpleExpression<String> civicSummaryHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return data.getValue(SUMMARY_URL.getName());
            }
        };
    }
}
