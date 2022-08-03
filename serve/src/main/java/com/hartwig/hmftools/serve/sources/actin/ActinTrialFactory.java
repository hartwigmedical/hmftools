package com.hartwig.hmftools.serve.sources.actin;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.cancertype.CancerTypeConstants;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.treatment.ImmutableTreatment;
import com.hartwig.hmftools.serve.treatment.Treatment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActinTrialFactory.class);

    private static final String TRIAL_COHORT_SEPARATOR = "|";

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry entry, @NotNull String sourceEvent) {
        return ImmutableActinTrial.builder()
                .source(Knowledgebase.ACTIN)
                .sourceEvent(sourceEvent)
                .sourceUrls(Sets.newHashSet())
                .treatment(extractTreatment(entry))
                .applicableCancerType(CancerTypeConstants.CANCER_TYPE)
                .blacklistCancerTypes(Sets.newHashSet())
                .level(EvidenceLevel.B)
                .direction(entry.isUsedAsInclusion() ? EvidenceDirection.RESPONSIVE : EvidenceDirection.NO_BENEFIT)
                .evidenceUrls(Sets.newHashSet())
                .build();
    }

    @NotNull
    private static Treatment extractTreatment(@NotNull ActinEntry entry) {
        String addon = Strings.EMPTY;
        if (entry.cohort() != null) {
            if (entry.cohort().contains(TRIAL_COHORT_SEPARATOR)) {
                LOGGER.warn("ACTIN entry cohort contains cohort separator: {}", entry);
            }
            addon = TRIAL_COHORT_SEPARATOR + entry.cohort();
        }

        return ImmutableTreatment.builder()
                .treament(entry.trial() + addon)
                .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                .relevantTreatmentApproaches(Sets.newHashSet())
                .build();
    }
}