package com.hartwig.hmftools.common.actionability.cnv;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CopyNumberEvidenceAnalyzer {
    private static final double REL_GAIN = 3.0;
    private static final double ABS_LOSS = 0.5;

    @NotNull
    private final List<ActionableCopyNumber> actionableCopyNumbers;

    CopyNumberEvidenceAnalyzer(@NotNull final List<ActionableCopyNumber> actionableCopyNumbers) {
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
    public List<EvidenceItem> evidenceForCopyNumberEvent(@NotNull GeneCopyNumber geneCopyNumber, @Nullable String doidsPrimaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer, final double purplePloidy) {
        Double minCopyValue = (double) Math.max(0, Math.round(geneCopyNumber.minCopyNumber()));
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (ActionableCopyNumber actionableCopyNumber : actionableCopyNumbers) {
            if (checkCNVType(minCopyValue, purplePloidy).equals(actionableCopyNumber.cnvType())) {
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
    private static String checkCNVType(final double copyNumber, final double purplePloidy) {
        Double relativeCopyNumber = copyNumber / purplePloidy;
        String CNVType = "";
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            CNVType = "Deletion";
        } else if (Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN)) {
            CNVType = "Amplification";
        }
        return CNVType;
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableCopyNumber(@NotNull ActionableCopyNumber actionableCopyNumber) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableCopyNumber.reference())
                .source(actionableCopyNumber.source())
                .drug(actionableCopyNumber.drug())
                .drugsType(actionableCopyNumber.drugsType())
                .level(actionableCopyNumber.level())
                .response(actionableCopyNumber.response());
    }
}
