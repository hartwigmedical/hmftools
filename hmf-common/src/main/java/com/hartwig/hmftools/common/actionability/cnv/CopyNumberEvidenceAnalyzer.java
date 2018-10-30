package com.hartwig.hmftools.common.actionability.cnv;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

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
        Set<String> genes = Sets.newHashSet();
        for (ActionableCopyNumber cnvs : actionableCopyNumbers) {
            genes.add(cnvs.gene());
        }
        return genes;
    }

    @NotNull
    public List<EvidenceItem> evidenceForCopyNumber(@NotNull GeneCopyNumber geneCopyNumber, @Nullable String doidsPrimaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        // KODU: Assume the gene copy number has already been determined to be a significant event (LOSS or GAIN)
        for (ActionableCopyNumber actionableCopyNumber : actionableCopyNumbers) {
            if (checkCNVType(geneCopyNumber.value()).equals(actionableCopyNumber.cnvType()) && actionableCopyNumber.gene()
                    .equals(geneCopyNumber.gene())) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableCopyNumber(actionableCopyNumber);

                evidenceBuilder.event(geneCopyNumber.gene() + " " + actionableCopyNumber.cnvType());
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.foundTumorLocation(actionableCopyNumber.cancerType(),
                        doidsPrimaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }
        return evidenceItems;
    }

    @NotNull
    private static String checkCNVType(int copyNumber) {
        return copyNumber <= 1 ? "Deletion" : "Amplification";
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableCopyNumber(@NotNull ActionableCopyNumber actionableCopyNumber) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableCopyNumber.reference())
                .source(ActionabilitySource.fromString(actionableCopyNumber.source()))
                .drug(actionableCopyNumber.drug())
                .drugsType(actionableCopyNumber.drugsType())
                .level(actionableCopyNumber.level())
                .response(actionableCopyNumber.response());
    }
}
