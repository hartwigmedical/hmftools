package com.hartwig.hmftools.protect.curation;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;

import org.jetbrains.annotations.NotNull;

public final class EvidenceCuration {

    private static final Set<String> EVENT_REPORTING_BLACKLIST_KEYWORDS = Sets.newHashSet("TP53");
    private static final Set<String> TREATMENT_REPORTING_BLACKLIST = Sets.newHashSet("Chemotherapy", "Aspirin", "Steroids");
    @VisibleForTesting
    static final Set<CurationKey> BLACKLIST_KEYS = Sets.newHashSet();

    static {
        BLACKLIST_KEYS.add(ImmutableCurationKey.builder()
                .eventKeyword("PIK3CA")
                .treatment("Trastuzumab")
                .direction(EvidenceDirection.RESPONSIVE)
                .build());
    }

    private EvidenceCuration() {
    }

    @NotNull
    public static List<ProtectEvidence> applyReportingBlacklist(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (hasBlacklistedEvent(evidence) || hasBlacklistedTreatment(evidence) || keyIsBlacklisted(evidence)) {
                result.add(ImmutableProtectEvidence.builder().from(evidence).reported(false).build());
            } else {
                result.add(evidence);
            }
        }
        return result;
    }

    private static boolean hasBlacklistedEvent(@NotNull ProtectEvidence evidence) {
        for (String entry : EVENT_REPORTING_BLACKLIST_KEYWORDS) {
            if (evidence.genomicEvent().contains(entry)) {
                return true;
            }
        }

        return false;
    }

    private static boolean hasBlacklistedTreatment(@NotNull ProtectEvidence evidence) {
        for (String entry : TREATMENT_REPORTING_BLACKLIST) {
            if (evidence.treatment().equals(entry)) {
                return true;
            }
        }
        return false;
    }

    private static boolean keyIsBlacklisted(@NotNull ProtectEvidence evidence) {
        for (CurationKey key : BLACKLIST_KEYS) {
            if (evidence.genomicEvent().contains(key.eventKeyword()) && evidence.treatment().equals(key.treatment()) && evidence.direction() == key.direction()) {
                return true;
            }
        }

        return false;
    }
}
