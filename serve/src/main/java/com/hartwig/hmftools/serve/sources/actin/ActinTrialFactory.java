package com.hartwig.hmftools.serve.sources.actin;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.tumorlocation.ImmutableTumorLocation;
import com.hartwig.hmftools.serve.tumorlocation.TumorLocation;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry actionTrial, @NotNull String rawInput) {

        Set<TumorLocation> tumorLocationBlacklistings = Sets.newHashSet();
        tumorLocationBlacklistings.add(ImmutableTumorLocation.builder()
                .cancerType("Hematologic cancer")
                .doid("2531")
                .build());

        return ImmutableActinTrial.builder()
                .source(Knowledgebase.ACTIN)
                .sourceEvent(rawInput)
                .sourceUrls(Sets.newHashSet())
                .treatment(actionTrial.trial())
                .whiteList(ImmutableTumorLocation.builder()
                        .cancerType("Advanced Solid Tumor")
                        .doid("162")
                        .build())
                .blacklistings(tumorLocationBlacklistings)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Sets.newHashSet())
                .build();
    }
}