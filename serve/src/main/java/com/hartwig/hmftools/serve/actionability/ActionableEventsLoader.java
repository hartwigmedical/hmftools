package com.hartwig.hmftools.serve.actionability;

import java.io.IOException;
import java.util.List;

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
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActionableEventsLoader {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEventsLoader.class);

    private ActionableEventsLoader() {
    }

    @NotNull
    public static ActionableEvents readFromDir(@NotNull String actionabilityDir, @NotNull RefGenomeVersion refGenomeVersion)
            throws IOException {
        LOGGER.info("Loading SERVE actionability files from {}", actionabilityDir);

        String actionableHotspotTsv = ActionableHotspotFile.actionableHotspotTsvPath(actionabilityDir, refGenomeVersion);
        List<ActionableHotspot> hotspots = ActionableHotspotFile.read(actionableHotspotTsv);
        LOGGER.info(" Loaded {} actionable hotspots from {}", hotspots.size(), actionableHotspotTsv);

        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvPath(actionabilityDir, refGenomeVersion);
        List<ActionableRange> ranges = ActionableRangeFile.read(actionableRangeTsv);
        LOGGER.info(" Loaded {} actionable ranges from {}", ranges.size(), actionableRangeTsv);

        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvPath(actionabilityDir, refGenomeVersion);
        List<ActionableGene> genes = ActionableGeneFile.read(actionableGeneTsv);
        LOGGER.info(" Loaded {} actionable genes from {}", genes.size(), actionableGeneTsv);

        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(actionabilityDir, refGenomeVersion);
        List<ActionableFusion> fusions = ActionableFusionFile.read(actionableFusionTsv);
        LOGGER.info(" Loaded {} actionable fusions from {}", fusions.size(), actionableFusionTsv);

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(actionabilityDir, refGenomeVersion);
        List<ActionableSignature> signatures = ActionableSignatureFile.read(actionableSignatureTsv);
        LOGGER.info(" Loaded {} actionable signatures from {}", signatures.size(), actionableSignatureTsv);

        return ImmutableActionableEvents.builder()
                .hotspots(hotspots)
                .ranges(ranges)
                .genes(genes)
                .fusions(fusions)
                .signatures(signatures)
                .build();
    }
}
