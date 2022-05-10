package com.hartwig.hmftools.common.protect;

import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ProtectTestFactory {

    private ProtectTestFactory() {
    }

    @NotNull
    public static ImmutableProtectEvidence.Builder builder() {
        Set<ProtectSource> source = Sets.newHashSet();
        source.add(ImmutableProtectSource.builder()
                .name(Knowledgebase.CKB)
                .sourceEvent("hotspot")
                .sourceUrls(Sets.newHashSet())
                .evidenceType(ProtectEvidenceType.ANY_MUTATION)
                        .evidenceUrls(Sets.newHashSet())
                .build());

        return ImmutableProtectEvidence.builder()
                .event(Strings.EMPTY)
                .germline(false)
                .reported(true)
                .treatment(Strings.EMPTY)
                .onLabel(false)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(source);
    }
}