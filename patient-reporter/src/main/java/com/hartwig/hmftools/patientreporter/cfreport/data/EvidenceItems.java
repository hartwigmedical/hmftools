package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class EvidenceItems {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceItems.class);

    private EvidenceItems() {
    }

    @NotNull
    public static List<ProtectEvidence> sort(@NotNull List<ProtectEvidence> evidenceItems) {
        return evidenceItems.stream().sorted((item1, item2) -> {
            if (item1.level().equals(item2.level())) {
                if (item1.genomicEvent().equals(item2.genomicEvent())) {
                    return item1.genomicEvent().compareTo(item2.genomicEvent());
                } else {
                    return item1.genomicEvent().compareTo(item2.genomicEvent());
                }
            } else {
                return item1.level().name().compareTo(item2.level().name());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String sources(@NotNull ProtectEvidence evidence) {
        StringJoiner sourceString = new StringJoiner(",");
        for (Knowledgebase source: evidence.sources()) {
            sourceString.add(source.reportDisplay());
        }

        return sourceString.toString();
    }

    @NotNull
    public static String evidenceUrl(@NotNull ProtectEvidence evidence) {
        if (evidence.urls().isEmpty()) {
            LOGGER.warn("No URL configured for evidence '{}'", evidence);
            return Strings.EMPTY;
        }

        // We prefer pubmed URLs over all other URLs so if there is one pubmed then we use that.
        for (String url : evidence.urls()) {
            if (url.contains("pubmed")) {
                return url;
            }
        }

        // If there are no pubmeds, and the first url refers to google we remove it.
        String url = evidence.urls().iterator().next();
        if (url.contains("google")) {
            return Strings.EMPTY;
        } else {
            return url;
        }
    }

    public static int uniqueEventCount(@NotNull List<ProtectEvidence> evidenceItems) {
        Set<String> events = Sets.newHashSet();
        for (ProtectEvidence evidence : evidenceItems) {
            events.add(evidence.genomicEvent());
        }
        return events.size();
    }

    public static int uniqueTherapyCount(@NotNull List<ProtectEvidence> evidenceItems) {
        Set<String> treatments = Sets.newHashSet();
        for (ProtectEvidence evidence : evidenceItems) {
            treatments.add(evidence.treatment());
        }
        return treatments.size();
    }
}
