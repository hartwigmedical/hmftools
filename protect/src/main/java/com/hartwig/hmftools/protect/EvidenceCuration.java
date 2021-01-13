package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;

public final class EvidenceCuration {

    private static final Set<String> EVENT_REPORTING_BLACKLIST_KEYWORDS = Sets.newHashSet("TP53");
    private static final Set<String> TREATMENT_REPORTING_BLACKLIST = Sets.newHashSet("Chemotherapy", "Aspirin");

    private EvidenceCuration() {
    }

    @NotNull
    public static List<ProtectEvidence> applyReportingBlacklist(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (hasBlacklistedEvent(evidence) || hasBlacklistedTreatment(evidence)) {
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
}
