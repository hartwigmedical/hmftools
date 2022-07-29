package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RelevantTreatmentAprroachCuration {

    private static final Logger LOGGER = LogManager.getLogger(RelevantTreatmentAprroachCuration.class);

    @NotNull
    private final Map<RelevantTreatmentApprochCurationEntryKey, RelevantTreatmentApprochCurationEntry> curations;
    @NotNull
    private final Set<RelevantTreatmentApprochCurationEntry> usedCurations = Sets.newHashSet();

    public RelevantTreatmentAprroachCuration(
            @NotNull final Map<RelevantTreatmentApprochCurationEntryKey, RelevantTreatmentApprochCurationEntry> curations) {
        this.curations = curations;
    }

    public void reportUnusedFilterEntries() {
        int unusedFilterEntryCount = 0;
        for (Map.Entry<RelevantTreatmentApprochCurationEntryKey, RelevantTreatmentApprochCurationEntry> entry : curations.entrySet()) {
            LOGGER.debug("entry debug: " + entry);
            if (!usedCurations.contains(entry.getValue().curationKey())) {
                unusedFilterEntryCount++;
                LOGGER.warn(" Curation entry '{}' hasn't been used for treatment Approch curation", entry);
            }
        }

        LOGGER.debug(" Found {} unused filter entries during treatment Approch curation", unusedFilterEntryCount);
    }

    @NotNull
    public String isMatch(@NotNull RelevantTreatmentApprochCurationEntryKey key) {
        RelevantTreatmentApprochCurationEntry curationEntry = curations.get(key);

        if (curationEntry == null) {
            LOGGER.warn("The treatment '{}' with relevant treatment approach '{}' of event '{}' "
                            + "with level '{}' and direction '{}' isn't curated because missing in curation resource",
                    key.treatment(),
                    key.treatmentApproach(),
                    key.event(),
                    key.level(),
                    key.direction());
            return Strings.EMPTY;
        } else {
            switch (curationEntry.curationType()) {
                case EVENT_TREATMENT_APPROACH_CURATION_IGNORE: {
                    usedCurations.add(ImmutableRelevantTreatmentApprochCurationEntry.builder()
                            .curationType(RelevantTreatmentApproachCurationType.EVENT_TREATMENT_APPROACH_CURATION_IGNORE)
                            .curationKey(key)
                            .curatedtreatmentApproach(curationEntry.curatedtreatmentApproach())
                            .build());
                    return curationEntry.curatedtreatmentApproach() == null ? Strings.EMPTY : curationEntry.curatedtreatmentApproach();
                }
                case DIRECTION_TREATMENT_APPROACH_CURATION_IGNORE: {
                    usedCurations.add(ImmutableRelevantTreatmentApprochCurationEntry.builder()
                            .curationType(RelevantTreatmentApproachCurationType.DIRECTION_TREATMENT_APPROACH_CURATION_IGNORE)
                            .curationKey(key)
                            .curatedtreatmentApproach(curationEntry.curatedtreatmentApproach())
                            .build());
                    return curationEntry.curatedtreatmentApproach() == null ? Strings.EMPTY : curationEntry.curatedtreatmentApproach();
                }
                case TREATMENT_APPROACH_CURATION: {
                    usedCurations.add(ImmutableRelevantTreatmentApprochCurationEntry.builder()
                            .curationType(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION)
                            .curationKey(key)
                            .curatedtreatmentApproach(curationEntry.curatedtreatmentApproach())
                            .build());
                    return curationEntry.curatedtreatmentApproach() == null ? Strings.EMPTY : curationEntry.curatedtreatmentApproach();
                }
                default: {
                    LOGGER.warn("Curation entry found with unrecognized type: {}", curationEntry);
                    return Strings.EMPTY;
                }
            }
        }
    }
}