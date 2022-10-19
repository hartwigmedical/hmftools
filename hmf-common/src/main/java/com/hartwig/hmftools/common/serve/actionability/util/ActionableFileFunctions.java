package com.hartwig.hmftools.common.serve.actionability.util;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.serve.actionability.ImmutableTreatment;
import com.hartwig.hmftools.common.serve.actionability.Treatment;
import com.hartwig.hmftools.common.serve.cancertype.CancerType;
import com.hartwig.hmftools.common.serve.cancertype.CancerTypeFactory;
import com.hartwig.hmftools.common.serve.cancertype.ImmutableCancerType;

import org.jetbrains.annotations.NotNull;

public final class ActionableFileFunctions {

    public static final String FIELD_DELIMITER = "\t";

    private static final String URL_DELIMITER = ",";

    private ActionableFileFunctions() {
    }

    @NotNull
    public static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("source")
                .add("sourceEvent")
                .add("sourceUrls")
                .add("treatment")
                .add("sourceRelevantTreatmentApproaches")
                .add("relevantTreatmentApproaches")
                .add("applicableCancerType")
                .add("applicableDoid")
                .add("blacklistCancerTypes")
                .add("level")
                .add("direction")
                .add("evidenceUrls")
                .toString();
    }

    @NotNull
    public static ActionableEvent fromLine(@NotNull String[] values, int startingPosition) {

        return new ActionableEvent() {

            @NotNull
            @Override
            public Knowledgebase source() {
                return Knowledgebase.valueOf(values[startingPosition]);
            }

            @NotNull
            @Override
            public String sourceEvent() {
                return values[startingPosition + 1];
            }

            @NotNull
            @Override
            public Set<String> sourceUrls() {
                int urlPosition = startingPosition + 2;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }

            @NotNull
            @Override
            public Treatment treatment() {
                int sourceDrugClassPosition = startingPosition + 4;
                int drugClassPosition = startingPosition + 5;
                return ImmutableTreatment.builder()
                        .treament(values[startingPosition + 3])
                        .sourceRelevantTreatmentApproaches(
                                values.length > sourceDrugClassPosition ? stringToDrugClasses(values[sourceDrugClassPosition]) : Sets.newHashSet())
                        .relevantTreatmentApproaches(
                                values.length > drugClassPosition ? stringToDrugClasses(values[drugClassPosition]) : Sets.newHashSet())
                        .build();

            }

            @NotNull
            @Override
            public CancerType applicableCancerType() {
                return ImmutableCancerType.builder().name(values[startingPosition + 6]).doid(values[startingPosition + 7]).build();
            }

            @NotNull
            @Override
            public Set<CancerType> blacklistCancerTypes() {
                return CancerTypeFactory.fromString(values[startingPosition + 8]);
            }

            @NotNull
            @Override
            public EvidenceLevel level() {
                return EvidenceLevel.valueOf(values[startingPosition + 9]);
            }

            @NotNull
            @Override
            public EvidenceDirection direction() {
                return EvidenceDirection.valueOf(values[startingPosition + 10]);
            }

            @NotNull
            @Override
            public Set<String> evidenceUrls() {
                int urlPosition = startingPosition + 11;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }
        };
    }

    @NotNull
    public static String toLine(@NotNull ActionableEvent event) {
        return new StringJoiner(FIELD_DELIMITER).add(event.source().toString())
                .add(event.sourceEvent())
                .add(urlsToString(event.sourceUrls()))
                .add(event.treatment().treament())
                .add(drugClassesToString(event.treatment().sourceRelevantTreatmentApproaches()))
                .add(drugClassesToString(event.treatment().relevantTreatmentApproaches()))
                .add(event.applicableCancerType().name())
                .add(event.applicableCancerType().doid())
                .add(CancerTypeFactory.toString(event.blacklistCancerTypes()))
                .add(event.level().toString())
                .add(event.direction().toString())
                .add(urlsToString(event.evidenceUrls()))
                .toString();
    }

    @NotNull
    private static Set<String> stringToUrls(@NotNull String fieldValue) {
        return Sets.newHashSet(fieldValue.split(URL_DELIMITER));
    }

    @NotNull
    public static String urlsToString(@NotNull Set<String> urls) {
        StringJoiner joiner = new StringJoiner(URL_DELIMITER);
        for (String url : urls) {
            joiner.add(url);
        }
        return joiner.toString();
    }

    @NotNull
    private static Set<String> stringToDrugClasses(@NotNull String fieldValue) {
        return Sets.newHashSet(fieldValue.split(URL_DELIMITER));
    }

    @NotNull
    public static String drugClassesToString(@NotNull Set<String> drugClasses) {
        StringJoiner joiner = new StringJoiner(URL_DELIMITER);
        for (String url : drugClasses) {
            joiner.add(url);
        }
        return joiner.toString();
    }
}