package com.hartwig.hmftools.common.actionability.cnv;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CopyNumberEvidenceAnalyzer {

    @NotNull
    private final List<ActionableCopyNumber> actionableCopyNumbers;

    CopyNumberEvidenceAnalyzer(@NotNull List<ActionableCopyNumber> actionableCopyNumbers) {
        this.actionableCopyNumbers = actionableCopyNumbers;
    }

    @NotNull
    public Set<String> actionableGenes() {
        return actionableCopyNumbers.stream().map(ActionableCopyNumber::gene).collect(Collectors.toSet());
    }

    @NotNull
    public List<ActionableCopyNumber> actionableCopyNumbers() {
        return actionableCopyNumbers;
    }

    @NotNull
    public List<EvidenceItem> evidenceForCopyNumber(@NotNull ReportableGainLoss reportableGainLoss, @Nullable String primaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (ActionableCopyNumber actionableCopyNumber : actionableCopyNumbers) {
            if (typeMatches(reportableGainLoss, actionableCopyNumber) && actionableCopyNumber.gene().equals(reportableGainLoss.gene())) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableCopyNumber(actionableCopyNumber);
                evidenceBuilder.event(reportableGainLoss.gene() + " " + actionableCopyNumber.type().readableString());
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.isCancerTypeMatch(actionableCopyNumber.cancerType(), primaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }

        return evidenceItems;
    }

    private static boolean typeMatches(@NotNull ReportableGainLoss reportableGainLoss, @NotNull ActionableCopyNumber actionableCopyNumber) {
        CopyNumberType geneType = reportableGainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                || reportableGainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                ? CopyNumberType.DELETION
                : CopyNumberType.AMPLIFICATION;
        return geneType == actionableCopyNumber.type();
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableCopyNumber(@NotNull ActionableCopyNumber actionableCopyNumber) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableCopyNumber.reference())
                .source(ActionabilitySource.fromString(actionableCopyNumber.source()))
                .drug(actionableCopyNumber.drug())
                .drugsType(actionableCopyNumber.drugsType())
                .level(EvidenceLevel.fromString(actionableCopyNumber.level()))
                .response(actionableCopyNumber.response())
                .cancerType(actionableCopyNumber.cancerType())
                .scope(EvidenceScope.SPECIFIC);
    }
}
