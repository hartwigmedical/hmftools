package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;

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
                    .protectSources(Sets.newHashSet())
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
                    .protectSources(entry.getValue().protectSources())
                    .build());
        }
        return consolidatedEvents;
    }

    private static class ConsolidatedData {

        @NotNull
        private final Set<ProtectSource> protectSources = Sets.newHashSet();

        public ConsolidatedData() {
        }

        public void appendEvidence(@NotNull ProtectEvidence evidence) {
            protectSources.addAll(evidence.protectSources());
        }

        @NotNull
        public Set<ProtectSource> protectSources() {
            return protectSources;
        }
    }
}