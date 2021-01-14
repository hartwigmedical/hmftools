package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class EvidenceConsolidation {

    private EvidenceConsolidation() {
    }

    @NotNull
    public static List<ProtectEvidence> consolidate(@NotNull List<ProtectEvidence> evidences) {
        Map<ProtectEvidence, ConsolidatedData> dataPerEvidence = Maps.newHashMap();
        for (ProtectEvidence evidence : evidences) {
            ProtectEvidence strippedEvidence =
                    ImmutableProtectEvidence.builder().from(evidence).urls(Sets.newHashSet()).sources(Sets.newHashSet()).build();
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
                    .sources(entry.getValue().sources())
                    .urls(entry.getValue().urls())
                    .build());
        }
        return consolidatedEvents;
    }

    private static class ConsolidatedData {

        @NotNull
        private final Set<Knowledgebase> sources = Sets.newHashSet();
        @NotNull
        private final Set<String> urls = Sets.newHashSet();

        public ConsolidatedData() {
        }

        public void appendEvidence(@NotNull ProtectEvidence evidence) {
            sources.addAll(evidence.sources());
            urls.addAll(evidence.urls());
        }

        @NotNull
        public Set<Knowledgebase> sources() {
            return sources;
        }

        @NotNull
        public Set<String> urls() {
            return urls;
        }
    }
}
