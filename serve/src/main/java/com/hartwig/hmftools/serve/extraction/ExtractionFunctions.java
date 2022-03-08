package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristicUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLAUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.actionability.util.ActionableEventUrlMerger;
import com.hartwig.hmftools.serve.extraction.codon.CodonFunctions;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonFunctions;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ExtractionFunctions {

    private static final Logger LOGGER = LogManager.getLogger(ExtractionFunctions.class);

    private ExtractionFunctions() {
    }

    @NotNull
    public static ExtractionResult consolidateActionableEvents(@NotNull ExtractionResult result) {
        return ImmutableExtractionResult.builder()
                .from(result)
                .actionableHotspots(ActionableEventUrlMerger.merge(result.actionableHotspots(), new ActionableHotspotUrlConsolidator()))
                .actionableRanges(ActionableEventUrlMerger.merge(result.actionableRanges(), new ActionableRangeUrlConsolidator()))
                .actionableGenes(ActionableEventUrlMerger.merge(result.actionableGenes(), new ActionableGeneUrlConsolidator()))
                .actionableFusions(ActionableEventUrlMerger.merge(result.actionableFusions(), new ActionableFusionUrlConsolidator()))
                .actionableCharacteristics(ActionableEventUrlMerger.merge(result.actionableCharacteristics(),
                        new ActionableCharacteristicUrlConsolidator()))
                .actionableHLA(ActionableEventUrlMerger.merge(result.actionableHLA(), new ActionableHLAUrlConsolidator()))
                .build();
    }

    @NotNull
    public static ExtractionResult merge(@NotNull List<ExtractionResult> results) {
        RefGenomeVersion version = uniqueVersion(results);
        ImmutableExtractionResult.Builder mergedBuilder = ImmutableExtractionResult.builder().refGenomeVersion(version);

        Set<KnownHotspot> allHotspots = Sets.newHashSet();
        Set<KnownCodon> allCodons = Sets.newHashSet();
        Set<KnownExon> allExons = Sets.newHashSet();
        Set<KnownCopyNumber> allCopyNumbers = Sets.newHashSet();
        Set<KnownFusionPair> allFusionPairs = Sets.newHashSet();

        for (ExtractionResult result : results) {

            Set<ActionableRange> actionableRange = curate(result.actionableRanges());

            allHotspots.addAll(result.knownHotspots());
            allCodons.addAll(result.knownCodons());
            allExons.addAll(result.knownExons());
            allCopyNumbers.addAll(result.knownCopyNumbers());
            allFusionPairs.addAll(result.knownFusionPairs());

            mergedBuilder.addAllActionableHotspots(result.actionableHotspots());
            mergedBuilder.addAllActionableRanges(actionableRange);
            mergedBuilder.addAllActionableGenes(result.actionableGenes());
            mergedBuilder.addAllActionableFusions(result.actionableFusions());
            mergedBuilder.addAllActionableCharacteristics(result.actionableCharacteristics());
            mergedBuilder.addAllActionableHLA(result.actionableHLA());
        }

        ExtractionResult mergedResult = mergedBuilder.knownHotspots(HotspotFunctions.consolidate(allHotspots))
                .knownCodons(CodonFunctions.consolidate(allCodons))
                .knownExons(ExonFunctions.consolidate(allExons))
                .knownCopyNumbers(CopyNumberFunctions.consolidate(allCopyNumbers))
                .knownFusionPairs(FusionFunctions.consolidate(allFusionPairs))
                .build();

        return consolidateActionableEvents(mergedResult);
    }

    @NotNull
    static Set<ActionableRange> curate(@NotNull Set<ActionableRange> actionableRanges) {
        Set<ActionableRange> curatedRanges = Sets.newHashSet();
        for (ActionableRange range : actionableRanges) {
            if (range.gene().equals("BRAF") && range.rank() == 600 && range.rangeType().equals(RangeType.CODON)) {
                int start = 140753335;
                int end = 140753337;
                String canonicalTranscriptID = "ENST00000288602";
                curatedRanges.add(ImmutableActionableRange.builder()
                        .from(range)
                        .start(start)
                        .end(end)
                        .transcript(canonicalTranscriptID)
                        .build());
            } else {
                curatedRanges.add(range);
            }
        }
        return curatedRanges;
    }

    @NotNull
    private static RefGenomeVersion uniqueVersion(@NotNull List<ExtractionResult> results) {
        if (results.isEmpty()) {
            RefGenomeVersion defaultVersion = RefGenomeVersion.V38;
            LOGGER.warn("Cannot extract ref genome version for empty list of results. Reverting to default {}", defaultVersion);
            return defaultVersion;
        }

        RefGenomeVersion version = results.get(0).refGenomeVersion();
        for (ExtractionResult result : results) {
            if (result.refGenomeVersion() != version) {
                throw new IllegalStateException("Ref genome version is not unique amongst list of extraction results");
            }
        }

        return version;
    }
}