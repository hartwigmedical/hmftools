package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RelevantTreatmentAprroachCurationTest {

    @Test
    public void canTestMatchEntries() {
        Map<RelevantTreatmentApprochCurationEntryKey, RelevantTreatmentApprochCurationEntry> curationEntries = Maps.newHashMap();
        curationEntries.put(canGenerateCurationKey("A", "A", "BRAF amplification", EvidenceLevel.A, EvidenceDirection.RESPONSIVE),
                canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION,
                        "A",
                        "A",
                        "BRAF amplification",
                        EvidenceLevel.A,
                        EvidenceDirection.RESPONSIVE,
                        "AA"));

        curationEntries.put(canGenerateCurationKey("B", "B", "BRAF amplification", EvidenceLevel.A, EvidenceDirection.RESPONSIVE),
                canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION,
                        "B",
                        "B",
                        "BRAF amplification",
                        EvidenceLevel.A,
                        EvidenceDirection.RESPONSIVE,
                        "BB"));

        curationEntries.put(canGenerateCurationKey("C", "C", "BRAF amplification", EvidenceLevel.A, EvidenceDirection.RESPONSIVE),
                canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION_IGNORE,
                        "C",
                        "C",
                        "BRAF amplification",
                        EvidenceLevel.A,
                        EvidenceDirection.RESPONSIVE,
                        Strings.EMPTY));

        RelevantTreatmentAprroachCuration curator = new RelevantTreatmentAprroachCuration(curationEntries);

        RelevantTreatmentApprochCurationEntryKey keyMatch1 = ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                .treatment("A")
                .treatmentApproach("A")
                .event("BRAF amplification")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();

        RelevantTreatmentApprochCurationEntryKey keyIgnore = ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                .treatment("C")
                .treatmentApproach("C")
                .event("BRAF amplification")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();

        RelevantTreatmentApprochCurationEntryKey keyUnmatch = ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                .treatment("D")
                .treatmentApproach("D")
                .event("BRAF amplification")
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .build();

        assertEquals("AA", curator.isMatch(keyMatch1));
        assertEquals(Strings.EMPTY, curator.isMatch(keyIgnore));
        assertEquals(Strings.EMPTY, curator.isMatch(keyUnmatch));
    }

    @NotNull
    public static RelevantTreatmentApprochCurationEntry canGenerateCurationEntry(@NotNull RelevantTreatmentApproachCurationType type,
            @NotNull String treatment, @NotNull String treatmentApproach, @NotNull String event, @NotNull EvidenceLevel level,
            @NotNull EvidenceDirection direction, @NotNull String curation) {
        return ImmutableRelevantTreatmentApprochCurationEntry.builder()
                .curationType(type)
                .curationKey(canGenerateCurationKey(treatment, treatmentApproach, event, level, direction))
                .curatedtreatmentApproach(curation)
                .build();
    }

    @NotNull
    public static RelevantTreatmentApprochCurationEntryKey canGenerateCurationKey(@NotNull String treatment,
            @NotNull String treatmentApproach, @NotNull String event, @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction) {
        return ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                .treatment(treatment)
                .treatmentApproach(treatmentApproach)
                .event(event)
                .level(level)
                .direction(direction)
                .build();
    }
}