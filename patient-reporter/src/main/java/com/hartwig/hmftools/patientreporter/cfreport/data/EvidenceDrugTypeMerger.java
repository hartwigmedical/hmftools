package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceItemMerger;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItemMerger;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class EvidenceDrugTypeMerger {

    private EvidenceDrugTypeMerger() {
    }

    @NotNull
    public static List<EvidenceItemMerger> merge(List<EvidenceItem> items) {
        StringBuilder drugsString = new StringBuilder();

        Map<DrugsKey, EvidenceItem> mapEvidence = Maps.newHashMap();
        for (EvidenceItem item : items) {
            mapEvidence.put(new DrugsKey(item.event(),
                    item.scope(),
                    item.level(),
                    item.response(),
                    item.drugsType(),
                    item.source()), item);
        }

        List<EvidenceItemMerger> evidenceItems = Lists.newArrayList();
        List<String> drug = Lists.newArrayList();
        for (Map.Entry<DrugsKey, EvidenceItem> entry : mapEvidence.entrySet()) {
            EvidenceItem itemsForKey = entry.getValue();
            if (mapEvidence.containsKey(entry.getKey())) {
                drug.add(itemsForKey.drug());
            }

            for (String drugs : drug) {
                drugsString.append(drugs).append(",");
            }

            evidenceItems.add(ImmutableEvidenceItemMerger.builder()
                    .event(itemsForKey.event())
                    .source(itemsForKey.source())
                    .reference(itemsForKey.reference())
                    .drug(itemsForKey.drugsType() + " (" + drugsString + ")")
                    .level(itemsForKey.level())
                    .response(itemsForKey.response())
                    .isOnLabel(itemsForKey.isOnLabel())
                    .cancerType(itemsForKey.cancerType())
                    .scope(itemsForKey.scope())
                    .build());
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
            return Objects.equals(event, drugsKey.event) && Objects.equals(match, drugsKey.match) && Objects.equals(level, drugsKey.level)
                    && Objects.equals(response, drugsKey.response) && Objects.equals(drugType, drugsKey.drugType) && Objects.equals(source,
                    drugsKey.source);
        }

        @Override
        public int hashCode() {
            return Objects.hash(event, match, level, response, drugType, source);
        }
    }
}
