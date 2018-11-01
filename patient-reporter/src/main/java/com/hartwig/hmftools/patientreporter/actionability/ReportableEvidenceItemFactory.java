package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableEvidenceItemFactory {

    private ReportableEvidenceItemFactory() {
    }

    @NotNull
    public static List<EvidenceItem> extractNonTrials(@NotNull List<EvidenceItem> evidenceItems) {
        return evidenceItems.stream().filter(evidenceItem -> !evidenceItem.source().isTrialSource()).collect(Collectors.toList());
    }

    @NotNull
    public static List<EvidenceItem> reportableFlatList(@NotNull Map<?, List<EvidenceItem>> evidenceItemMap) {
        return filterForReporting(toList(evidenceItemMap));
    }

    @NotNull
    private static List<EvidenceItem> toList(@NotNull Map<?, List<EvidenceItem>> evidenceItemMap) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> items : evidenceItemMap.values()) {
            evidenceItemList.addAll(items);
        }
        return evidenceItemList;
    }

    @NotNull
    private static List<EvidenceItem> filterForReporting(@NotNull List<EvidenceItem> evidenceItems) {
        Map<Key, List<EvidenceItem>> itemsPerKey = Maps.newHashMap();
        for (EvidenceItem item : evidenceItems) {
            Key key = new Key(item.event(), item.drug(), item.response());
            List<EvidenceItem> items = itemsPerKey.get(key);
            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(item);
            itemsPerKey.put(key, items);
        }

        List<EvidenceItem> evidenceFiltered = Lists.newArrayList();
        for (Map.Entry<Key, List<EvidenceItem>> entry : itemsPerKey.entrySet()) {
            List<EvidenceItem> itemsForKey = entry.getValue();
            EvidenceItem highestOnLabel = highestOnLabel(itemsForKey);
            EvidenceItem highestOffLabel = highestOffLabel(itemsForKey);

            if (highestOnLabel != null && highestOffLabel != null) {
                evidenceFiltered.add(highestOnLabel);
                if (hasHigherEvidence(highestOffLabel, highestOnLabel)) {
                    evidenceFiltered.add(highestOffLabel);
                }
            } else if (highestOnLabel != null) {
                evidenceFiltered.add(highestOnLabel);
            } else {
                evidenceFiltered.add(highestOffLabel);
            }
        }

        List<EvidenceItem> evidenceFilteredOnLevel = Lists.newArrayList();
        for (EvidenceItem evidence : evidenceFiltered) {
            if (hasReportableEvidenceLevel(evidence)) {
                evidenceFilteredOnLevel.add(evidence);
            }
        }

        return evidenceFilteredOnLevel;
    }

    @VisibleForTesting
    static boolean hasReportableEvidenceLevel(@NotNull EvidenceItem evidence) {
        return evidence.level().isReportedEvidenceItemLevel();
    }

    @VisibleForTesting
    static boolean hasHigherEvidence(@NotNull EvidenceItem item1, @NotNull EvidenceItem item2) {
        return item1.level().compareTo(item2.level()) < 0;
    }

    @VisibleForTesting
    @Nullable
    static EvidenceItem highestOffLabel(final List<EvidenceItem> items) {
        return highest(items, false);
    }

    @VisibleForTesting
    @Nullable
    static EvidenceItem highestOnLabel(@NotNull List<EvidenceItem> items) {
        return highest(items, true);
    }

    @Nullable
    private static EvidenceItem highest(@NotNull List<EvidenceItem> items, boolean shouldBeOnLabel) {
        EvidenceItem highest = null;
        for (EvidenceItem item : items) {
            if (item.isOnLabel() == shouldBeOnLabel) {
                if (highest == null || hasHigherEvidence(item, highest)) {
                    highest = item;
                }
            }
        }
        return highest;
    }

    private static class Key {

        @NotNull
        private final String event;
        @NotNull
        private final String drug;
        @NotNull
        private final String response;

        private Key(@NotNull final String event, @NotNull final String drug, @NotNull final String response) {
            this.event = event;
            this.drug = drug;
            this.response = response;
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
            return Objects.equals(event, key.event) && Objects.equals(drug, key.drug) && Objects.equals(response, key.response);
        }

        @Override
        public int hashCode() {
            return Objects.hash(event, drug, response);
        }
    }
}