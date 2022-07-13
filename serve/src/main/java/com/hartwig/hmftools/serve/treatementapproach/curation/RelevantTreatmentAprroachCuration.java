package com.hartwig.hmftools.serve.treatementapproach.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RelevantTreatmentAprroachCuration {

    private static final Logger LOGGER = LogManager.getLogger(RelevantTreatmentAprroachCuration.class);

    @NotNull
    private final List<RelevantTreatmentApprochCurationEntry> curations;
    @NotNull
    private final Set<RelevantTreatmentApprochCurationEntry> usedCurations = Sets.newHashSet();

    public RelevantTreatmentAprroachCuration(@NotNull final List<RelevantTreatmentApprochCurationEntry> curations) {
        this.curations = curations;
    }

    public void reportUnusedFilterEntries() {
        int unusedFilterEntryCount = 0;
        for (RelevantTreatmentApprochCurationEntry entry : curations) {
            if (!usedCurations.contains(entry.curationKey())) {
                unusedFilterEntryCount++;
                LOGGER.warn(" Curation entry '{}' hasn't been used for treatment Approch curation", entry);
            }
        }

        LOGGER.debug(" Found {} unused filter entries during treatment Approch curation", unusedFilterEntryCount);
    }

    @NotNull
    public String isMatch(
            @NotNull RelevantTreatmentApprochCurationEntryKey key) {
        for (RelevantTreatmentApprochCurationEntry curationEntry : curations) {
            switch (curationEntry.curationType()) {
                case TREATMENT_APPROACH_CURATION_IGNORE: {
                    if (curationEntry.curationKey() == key) {
                        usedCurations.add(ImmutableRelevantTreatmentApprochCurationEntry.builder()
                                .curationType(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION_IGNORE)
                                .curationKey(key)
                                .curatedtreatmentApproach(curationEntry.curatedtreatmentApproach())
                                .build());
                        return Strings.EMPTY;
                    }
                }
                case TREATMENT_APPROACH_CURATION: {
                    if (curationEntry.curationKey().equals(key)) {
                        usedCurations.add(ImmutableRelevantTreatmentApprochCurationEntry.builder()
                                .curationType(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION)
                                .curationKey(key)
                                .curatedtreatmentApproach(curationEntry.curatedtreatmentApproach())
                                .build());
                        return curationEntry.curatedtreatmentApproach();
                    } else {
                        LOGGER.warn("The treatment '{}' with relevant treatment approach '{}' of event '{}' "
                                        + "with level '{}' and direction '{}' isn't curated",
                                key.treatment(),
                                key.treatmentApproach(),
                                key.event(),
                                key.level(),
                                key.direction());
                        return Strings.EMPTY;
                    }
                }
                default: {
                    LOGGER.warn("Curation entry found with unrecognized type: {}", curationEntry);
                    return Strings.EMPTY;
                }
            }
        }
        return Strings.EMPTY;
    }
}