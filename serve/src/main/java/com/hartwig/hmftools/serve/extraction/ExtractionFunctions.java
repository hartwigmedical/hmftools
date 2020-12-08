package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;

import org.jetbrains.annotations.NotNull;

public final class ExtractionFunctions {

    private ExtractionFunctions() {
    }

    @NotNull
    public static ExtractionResult merge(@NotNull List<ExtractionResult> results) {
        ImmutableExtractionResult.Builder mergedBuilder = ImmutableExtractionResult.builder();

        Set<KnownHotspot> allHotspots = Sets.newHashSet();
        Set<KnownCopyNumber> allCopyNumbers = Sets.newHashSet();
        Set<KnownFusionPair> allFusionPairs = Sets.newHashSet();

        for (ExtractionResult result : results) {
            allHotspots.addAll(result.knownHotspots());
            allCopyNumbers.addAll(result.knownCopyNumbers());
            allFusionPairs.addAll(result.knownFusionPairs());

            mergedBuilder.addAllActionableHotspots(result.actionableHotspots());
            mergedBuilder.addAllActionableRanges(result.actionableRanges());
            mergedBuilder.addAllActionableGenes(result.actionableGenes());
            mergedBuilder.addAllActionableFusions(result.actionableFusions());
            mergedBuilder.addAllActionableSignatures(result.actionableSignatures());
        }

        return mergedBuilder.knownHotspots(HotspotFunctions.consolidate(allHotspots))
                .knownCopyNumbers(CopyNumberFunctions.consolidate(allCopyNumbers))
                .knownFusionPairs(FusionFunctions.consolidate(allFusionPairs))
                .build();
    }
}
