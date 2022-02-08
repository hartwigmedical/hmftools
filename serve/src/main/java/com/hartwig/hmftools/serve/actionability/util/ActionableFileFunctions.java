package com.hartwig.hmftools.serve.actionability.util;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionableFileFunctions {

    public static final String FIELD_DELIMITER = "\t";

    private static final String URL_DELIMITER = ",";

    private ActionableFileFunctions() {
    }

    @NotNull
    public static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("rawInput")
                .add("source")
                .add("treatment")
                .add("cancerType")
                .add("doid")
                .add("blacklistCancerType")
                .add("blacklistedDoid")
                .add("level")
                .add("direction")
                .add("urlSource")
                .add("urls")
                .toString();
    }

    @NotNull
    public static ActionableEvent fromLine(@NotNull String[] values, int startingPosition) {
        return new ActionableEvent() {
            @NotNull
            @Override
            public String rawInput() {
                return values[startingPosition];
            }

            @NotNull
            @Override
            public Knowledgebase source() {
                return Knowledgebase.valueOf(values[startingPosition + 1]);
            }

            @NotNull
            @Override
            public String treatment() {
                return values[startingPosition + 2];
            }

            @NotNull
            @Override
            public String cancerType() {
                return values[startingPosition + 3];
            }

            @NotNull
            @Override
            public String doid() {
                return values[startingPosition + 4];
            }

            @NotNull
            @Override
            public Set<String> blacklistCancerType() {
                int urlPosition = startingPosition + 5;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }

            @NotNull
            @Override
            public Set<String> blacklistedDoid() {
                int urlPosition = startingPosition + 6;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
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
            public Set<String> urlSource(){
                int urlPosition = startingPosition + 9;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }

            @NotNull
            @Override
            public Set<String> urls() {
                int urlPosition = startingPosition + 10;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }
        };
    }

    @NotNull
    public static String toLine(@NotNull ActionableEvent event) {
        return new StringJoiner(FIELD_DELIMITER).add(event.rawInput())
                .add(event.source().toString())
                .add(event.treatment())
                .add(event.cancerType())
                .add(event.doid())
                .add(urlsToString(event.blacklistCancerType()))
                .add(urlsToString(event.blacklistedDoid()))
                .add(event.level().toString())
                .add(event.direction().toString())
                .add(urlsToString(event.urlSource()))
                .add(urlsToString(event.urls()))
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
