package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class ReportableEvidenceItemFactory {

    private ReportableEvidenceItemFactory() {
    }

    @NotNull
    public static List<ProtectEvidence> extractNonTrialsOnLabel(@NotNull List<ProtectEvidence> evidenceItems) {
        List<ProtectEvidence> nonTrials = Lists.newArrayList();
        for (ProtectEvidence evidence: evidenceItems) {
            for (ProtectSource source: evidence.protectSources()) {
                if (source.sources() != Knowledgebase.ICLUSION && evidence.onLabel()){
                    nonTrials.add(evidence);
                }
            }
        }
        return nonTrials;
    }

    @NotNull
    public static List<ProtectEvidence> extractNonTrialsOffLabel(@NotNull List<ProtectEvidence> evidenceItems) {
        List<ProtectEvidence> nonTrials = Lists.newArrayList();
        for (ProtectEvidence evidence: evidenceItems) {
            for (ProtectSource source: evidence.protectSources()) {
                if (source.sources() != Knowledgebase.ICLUSION && !evidence.onLabel()){
                    nonTrials.add(evidence);
                }
            }
        }
        return nonTrials;
    }
}