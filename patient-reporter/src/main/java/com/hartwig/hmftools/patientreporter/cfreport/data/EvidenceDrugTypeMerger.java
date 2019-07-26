package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceItemMerger;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItemMerger;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class EvidenceDrugTypeMerger {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceDrugTypeMerger.class);

    private EvidenceDrugTypeMerger() {
    }

    @NotNull
    public static List<EvidenceItemMerger> merge(List<EvidenceItem> items) {
        StringBuilder drugsStringBuilderMultipleKeys = new StringBuilder();
        StringBuilder drugsStringBuilder = new StringBuilder();

        String drugsStringTotal = Strings.EMPTY;

        List<EvidenceItemMerger> evidenceItems = Lists.newArrayList();

        Map<DrugsKey, List<EvidenceItem>> mapEvidence = Maps.newHashMap();
        for (EvidenceItem item : items) {

            DrugsKey drugsKey = new DrugsKey(item.event(), item.scope(), item.level(), item.response(), item.drugsType(), item.source());

            if (mapEvidence.containsKey(drugsKey) && !item.drugsType().equals("Unknown") && !item.drugsType().equals(Strings.EMPTY)) {
                List<EvidenceItem> itemsForKey = mapEvidence.get(drugsKey);
                itemsForKey.add(item);
                mapEvidence.put(drugsKey, itemsForKey);
                drugsStringBuilder.append(",").append(item.drug());
            } else {
                mapEvidence.put(drugsKey, Lists.newArrayList(item));
            }
        }

        Set<DrugsKey> keys = mapEvidence.keySet();
        List<DrugsKey> drugsKeys = Lists.newArrayList();
        for (Map.Entry<DrugsKey, List<EvidenceItem>> entry : mapEvidence.entrySet()) {
            List<EvidenceItem> itemsForKey = entry.getValue();

            for (EvidenceItem itemValue: itemsForKey) {
                if (drugsStringBuilder.length() > 0) {
                    drugsStringTotal = drugsStringBuilder.substring(1, drugsStringBuilder.length());
                }
                if (!drugsStringBuilder.toString().contains(itemValue.drug())) {
                    drugsStringBuilderMultipleKeys.append(itemValue.drug());
                    drugsStringTotal = drugsStringBuilderMultipleKeys.toString();

                }
                if (keys.contains(entry.getKey()) && !drugsKeys.contains(entry.getKey())) {
                    drugsKeys.add(entry.getKey());
                    LOGGER.info(entry.getKey());
                    LOGGER.info(drugsStringTotal);

                    evidenceItems.add(ImmutableEvidenceItemMerger.builder()
                            .event(itemValue.event())
                            .source(itemValue.source())
                            .reference(itemValue.reference())
                            .drug(itemValue.drugsType() + " (" + drugsStringTotal + ")")
                            .level(itemValue.level())
                            .response(itemValue.response())
                            .isOnLabel(itemValue.isOnLabel())
                            .cancerType(itemValue.cancerType())
                            .scope(itemValue.scope())
                            .build());
                }
            }

        }
        return evidenceItems;
    }

    private static class DrugsKey {

        @NotNull
        private final String event;

        @NotNull
        private final EvidenceScope match;

        @NotNull
        private final EvidenceLevel level;

        @NotNull
        private final String response;

        @NotNull
        private final String drugType;

        @NotNull
        private final ActionabilitySource source;

        public DrugsKey(@NotNull final String event, @NotNull final EvidenceScope match, @NotNull final EvidenceLevel level,
                @NotNull final String response, @NotNull final String drugType, @NotNull final ActionabilitySource source) {
            this.event = event;
            this.match = match;
            this.level = level;
            this.response = response;
            this.drugType = drugType;
            this.source = source;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final DrugsKey drugsKey = (DrugsKey) o;
            if (drugsKey.drugType.equals("Unknown") || drugsKey.drugType.equals(Strings.EMPTY)) {
                return false;
            }
            return Objects.equals(event, drugsKey.event) && match == drugsKey.match && level == drugsKey.level && Objects.equals(response,
                    drugsKey.response) && Objects.equals(drugType, drugsKey.drugType) && source == drugsKey.source;
        }

        @Override
        public int hashCode() {
            return Objects.hash(event, match, level, response, drugType, source);
        }
    }
}
