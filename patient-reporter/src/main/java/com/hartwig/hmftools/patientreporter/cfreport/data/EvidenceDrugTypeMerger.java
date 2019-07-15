package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;

public final class EvidenceDrugTypeMerger {

    private EvidenceDrugTypeMerger() {
    }

    @NotNull
    public static List<EvidenceItem> merge(List<EvidenceItem> items) {

        Map<Key, List<EvidenceItem>> mapEvidence = Maps.newHashMap();
        for (EvidenceItem item : items) {
            mapEvidence.put(new Key(item.event()), Lists.newArrayList(item));
        }
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        return evidenceItems;
    }

    private static class Key {

        @NotNull
        private final String event;

        public Key(@NotNull final String event) {
            this.event = event;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final Key key = (Key) o;
            return Objects.equals(event, key.event);
        }

        @Override
        public int hashCode() {
            return Objects.hash(event);
        }
    }
}
