package com.hartwig.hmftools.serve.sources.actin;

import java.util.List;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.blacklisting.ImmutableTumorLocationBlacklisting;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklist;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklisting;
import com.hartwig.hmftools.serve.sources.ImmutableSources;
import com.hartwig.hmftools.serve.sources.Sources;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry actionTrial, @NotNull String rawInput) {

        List<TumorLocationBlacklisting> tumorLocationBlacklistings = Lists.newArrayList();
        tumorLocationBlacklistings.add(ImmutableTumorLocationBlacklisting.builder()
                .blacklistCancerType("Hematologic cancer")
                .blacklistedDoid("2531")
                .build());
        String tumorLocationBlacklisting = TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings);

        Sources sources = ImmutableSources.builder().sourceEvent(rawInput).source(Knowledgebase.ACTIN).build();

        return ImmutableActinTrial.builder()
                .source(sources)
                .treatment(actionTrial.trial())
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .cancerType("Advanced Solid Tumor")
                .doid("162")
                .tumorLocationBlacklisting(tumorLocationBlacklisting)
                .sourceUrls(Sets.newHashSet())
                .evidenceUrls(Sets.newHashSet())
                .build();
    }
}