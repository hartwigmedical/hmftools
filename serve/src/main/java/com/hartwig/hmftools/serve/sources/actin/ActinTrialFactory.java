package com.hartwig.hmftools.serve.sources.actin;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.blacklisting.ImmutableTumorLocationBlacklisting;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklist;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklisting;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry actionTrial, @NotNull String rawInput) {

        Set<TumorLocationBlacklisting> tumorLocationBlacklistings = Sets.newHashSet();
        tumorLocationBlacklistings.add(ImmutableTumorLocationBlacklisting.builder()
                .blacklistCancerType("Hematologic cancer")
                .blacklistedDoid("2531")
                .build());
        String tumorLocationBlacklist = TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings);
        String tumorLocationBlacklistDoid = TumorLocationBlacklist.extractTumorLocationDoid(tumorLocationBlacklistings);

        return ImmutableActinTrial.builder()
                .rawInput(rawInput)
                .source(Knowledgebase.ACTIN)
                .treatment(actionTrial.trial())
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .cancerType("Advanced Solid Tumor")
                .doid("162")
                .blacklistCancerType(tumorLocationBlacklist)
                .blacklistedDoid(tumorLocationBlacklistDoid)
                .urls(Sets.newHashSet())
                .urlSource(Sets.newHashSet())
                .build();
    }
}
