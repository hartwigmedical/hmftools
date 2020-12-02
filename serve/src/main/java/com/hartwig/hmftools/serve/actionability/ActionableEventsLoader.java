package com.hartwig.hmftools.serve.actionability;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActionableEventsLoader {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEventsLoader.class);

    @NotNull
    public static ActionableEvents readFromDir(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion)
            throws IOException {
        LOGGER.info("Loading SERVE actionability files from {}", serveActionabilityDir);

        String actionableHotspotTsv = ActionableHotspotFile.actionableHotspotTsvPath(serveActionabilityDir, refGenomeVersion);
        List<ActionableHotspot> hotspots = ActionableHotspotFile.read(actionableHotspotTsv);
        LOGGER.info(" Actionable hotspots: {}", hotspots.size());

        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvPath(serveActionabilityDir, refGenomeVersion);
        List<ActionableRange> ranges = ActionableRangeFile.read(actionableRangeTsv);
        LOGGER.info(" Actionable ranges: {}", ranges.size());

        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvPath(serveActionabilityDir, refGenomeVersion);
        List<ActionableGene> genes = ActionableGeneFile.read(actionableGeneTsv);
        LOGGER.info(" Actionable genes: {}", genes.size());

        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(serveActionabilityDir, refGenomeVersion);
        List<ActionableFusion> fusions = ActionableFusionFile.read(actionableFusionTsv);
        LOGGER.info(" Actionable fusions: {}", fusions.size());

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(serveActionabilityDir, refGenomeVersion);
        List<ActionableSignature> signatures = ActionableSignatureFile.read(actionableSignatureTsv);
        LOGGER.info(" Actionable signatures: {}", signatures.size());

        return ImmutableActionableEvents.builder()
                .addAllHotspots(hotspots)
                .addAllRanges(ranges)
                .addAllGenes(genes)
                .addAllFusions(fusions)
                .addAllSignatures(signatures)
                .build();
    }
}
