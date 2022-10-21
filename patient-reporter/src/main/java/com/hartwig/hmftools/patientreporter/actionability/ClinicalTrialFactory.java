package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.KnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.serve.datamodel.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrialFactory {

    private ClinicalTrialFactory() {
    }

    @NotNull
    public static List<ProtectEvidence> extractOnLabelTrials(@NotNull List<ProtectEvidence> evidenceItems) {
        List<ProtectEvidence> trials = Lists.newArrayList();
        for (ProtectEvidence evidence : evidenceItems) {
            Set<KnowledgebaseSource> protectSources = Sets.newHashSet();
            for (KnowledgebaseSource protectSource: evidence.sources()) {
                if (protectSource.name() == Knowledgebase.ICLUSION && evidence.onLabel()) {
                    protectSources.add(protectSource);
                }
            }

            if (protectSources.size() >= 1) {
                trials.add(ImmutableProtectEvidence.builder().from(evidence).sources(protectSources).build());
            }

        }
        return trials;
    }
}