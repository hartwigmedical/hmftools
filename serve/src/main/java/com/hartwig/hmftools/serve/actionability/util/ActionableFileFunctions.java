package com.hartwig.hmftools.serve.actionability.util;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.tumorlocation.ImmutableTumorLocation;
import com.hartwig.hmftools.serve.tumorlocation.TumorLocation;
import com.hartwig.hmftools.serve.tumorlocation.TumorLocationFactory;

import org.apache.logging.log4j.util.Strings;
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
                .add("cancerType")
                .add("doid")
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
            public String treatment() {
                return values[startingPosition + 3];
            }

            @NotNull
            @Override
            public TumorLocation whiteList() {
                return ImmutableTumorLocation.builder().cancerType(values[startingPosition + 4]).doid(values[startingPosition + 5]).build();
            }

            @NotNull
            @Override
            public Set<TumorLocation> blacklistings() {
                return values[startingPosition + 6].equals(Strings.EMPTY) ? TumorLocationFactory.readTumorLocationBlacklistingString(values[
                        startingPosition + 6]): Sets.newHashSet();
            }

            @NotNull
            @Override
            public EvidenceLevel level() {
                return EvidenceLevel.valueOf(values[startingPosition + 7]);
            }

            @NotNull
            @Override
            public EvidenceDirection direction() {
                return EvidenceDirection.valueOf(values[startingPosition + 8]);
            }

            @NotNull
            @Override
            public Set<String> evidenceUrls() {
                int urlPosition = startingPosition + 9;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }
        };
    }

    @NotNull
    public static String toLine(@NotNull ActionableEvent event) {
        return new StringJoiner(FIELD_DELIMITER).add(event.source().toString())
                .add(event.sourceEvent())
                .add(urlsToString(event.sourceUrls()))
                .add(event.treatment())
                .add(event.whiteList().cancerType())
                .add(event.whiteList().doid())
                .add(TumorLocationFactory.extractTumorLocationBlacklisting(event.blacklistings()))
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
    private static String urlsToString(@NotNull Set<String> urls) {
        StringJoiner joiner = new StringJoiner(URL_DELIMITER);
        for (String url : urls) {
            joiner.add(url);
        }
        return joiner.toString();
    }
}