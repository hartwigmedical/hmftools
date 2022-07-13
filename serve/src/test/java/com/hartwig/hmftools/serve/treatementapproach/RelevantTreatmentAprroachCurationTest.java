package com.hartwig.hmftools.serve.treatementapproach;

import java.util.List;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RelevantTreatmentAprroachCurationTest {

    @Test
    public void canTestMatchEntries() {
        List<RelevantTreatmentApprochCurationEntry> curationEntries = Lists.newArrayList();
        curationEntries.add(canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION,
                "A",
                "A",
                "BRAF amplification",
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                "AA"));

        curationEntries.add(canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION,
                "B",
                "B",
                "BRAF amplification",
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                "BB"));

        curationEntries.add(canGenerateCurationEntry(RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION_IGNORE,
                "C",
                "C",
                "BRAF amplification",
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE, Strings.EMPTY));

        RelevantTreatmentAprroachCuration curator =
                new RelevantTreatmentAprroachCuration(curationEntries);

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
                .curationKey(ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                        .treatment(treatment)
                        .treatmentApproach(treatmentApproach)
                        .event(event)
                        .level(level)
                        .direction(direction)
                        .build())
                .curatedtreatmentApproach(curation)
                .build();
    }
}