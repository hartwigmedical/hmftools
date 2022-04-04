package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ImmutableProtectSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class EvidenceConsolidation {

    private EvidenceConsolidation() {
    }

    @NotNull
    public static List<ProtectEvidence> consolidate(@NotNull List<ProtectEvidence> evidences) {
        Map<ProtectEvidence, ConsolidatedData> dataPerEvidence = Maps.newHashMap();
        for (ProtectEvidence evidence : evidences) {
            ProtectEvidence strippedEvidence = ImmutableProtectEvidence.builder()
                    .from(evidence)
                    .evidenceUrls(Sets.newHashSet())
                    .protectSources(ImmutableProtectSource.builder().build())
                    .build();
            ConsolidatedData data = dataPerEvidence.get(strippedEvidence);
            if (data == null) {
                data = new ConsolidatedData();
            }
            data.appendEvidence(evidence);
            dataPerEvidence.put(strippedEvidence, data);
        }

        List<ProtectEvidence> consolidatedEvents = Lists.newArrayList();
        for (Map.Entry<ProtectEvidence, ConsolidatedData> entry : dataPerEvidence.entrySet()) {
            consolidatedEvents.add(ImmutableProtectEvidence.builder()
                    .from(entry.getKey())
                    .evidenceUrls(entry.getValue().evidenceUrls())
                    .protectSources(ImmutableProtectSource.builder()
                            .sources(entry.getValue().sources())
                            .sourceEvent(entry.getValue().sourceEvent())
                            .sourceUrls(entry.getValue().sourceUrls())
                            .build())
                    .build());
        }
        return consolidatedEvents;
    }

    private static class ConsolidatedData {

        @NotNull
        private final Set<String> evidenceUrls = Sets.newTreeSet();
        @NotNull
        private final Set<Knowledgebase> sources = Sets.newTreeSet();
        @NotNull
        private final Set<String> sourceEvent = Sets.newTreeSet();
        @NotNull
        private final Set<String> sourceUrls = Sets.newTreeSet();

        public ConsolidatedData() {
        }

        public void appendEvidence(@NotNull ProtectEvidence evidence) {
            evidenceUrls.addAll(evidence.evidenceUrls());
            sources.addAll(evidence.protectSources().sources());
            sourceEvent.addAll(evidence.protectSources().sourceEvent());
            sourceUrls.addAll(evidence.protectSources().sourceUrls());
        }

        @NotNull
        public Set<String> evidenceUrls() {
            return evidenceUrls;
        }

        @NotNull
        public Set<Knowledgebase> sources() {
            return sources;
        }

        @NotNull
        public Set<String> sourceEvent() {
            return sourceEvent;
        }

        @NotNull
        public Set<String> sourceUrls() {
            return sourceUrls;
        }
    }
}
