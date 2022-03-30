package com.hartwig.hmftools.serve.sources.actin;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry entry, @NotNull String sourceEvent) {
        return ImmutableActinTrial.builder()
                .source(Knowledgebase.ACTIN)
                .sourceEvent(sourceEvent)
                .sourceUrls(Sets.newHashSet())
                .treatment(entry.trial())
                .applicableCancerType(ImmutableCancerType.builder().name("Cancer").doid("162").build())
                .blacklistCancerTypes(Sets.newHashSet())
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Sets.newHashSet())
                .build();
    }
}