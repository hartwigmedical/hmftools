package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EvidenceDrugTypeMerger {
    private static final Logger LOGGER = LogManager.getLogger(EvidenceDrugTypeMerger.class);

    private EvidenceDrugTypeMerger() {
    }

    @NotNull
    public static List<EvidenceItem> merge(List<EvidenceItem> items) {

        Map<DrugsKey, List<EvidenceItem>> mapEvidence = Maps.newHashMap();
        for (EvidenceItem item : items) {
            mapEvidence.put(new DrugsKey(item.event(), item.scope(), item.level(), item.response(), item.drugsType()), Lists.newArrayList(item));
        }


        for (Map.Entry<DrugsKey, List<EvidenceItem>> entry : mapEvidence.entrySet()) {
            List<EvidenceItem> itemsForKey = entry.getValue();


        }


        List<EvidenceItem> evidenceItems = Lists.newArrayList();
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

        public DrugsKey(@NotNull final String event, @NotNull final EvidenceScope match, @NotNull final EvidenceLevel level,
                @NotNull final String response, @NotNull final String drugType) {
            this.event = event;
            this.match = match;
            this.level = level;
            this.response = response;
            this.drugType = drugType;
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
            return Objects.equals(event, drugsKey.event);
        }

        @Override
        public int hashCode() {
            return Objects.hash(event, match, level, response, drugType);
        }
    }
}
