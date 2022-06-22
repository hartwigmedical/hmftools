package com.hartwig.hmftools.serve.extraction;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.ImmutableActionableEvents;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumberFile;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPairFile;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspotFile;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class KnownEventsLoader {

    private static final Logger LOGGER = LogManager.getLogger(KnownEventsLoader.class);

    private KnownEventsLoader() {
    }

    @NotNull
    public static KnownEvents readFromDir(@NotNull String actionabilityDir, @NotNull RefGenomeVersion refGenomeVersion)
            throws IOException {
        LOGGER.info("Loading SERVE known files from {} using ref genome version '{}'", actionabilityDir, refGenomeVersion);

        String knownHotspotTsv = KnownHotspotFile.knownHotspotVcfPath(actionabilityDir, refGenomeVersion);
        List<KnownHotspot> hotspots = KnownHotspotFile.read(knownHotspotTsv);
        LOGGER.info(" Loaded {} knonw hotspots from {}", hotspots.size(), knownHotspotTsv);

        String knownCodonTsv = KnownCodonFile.knownCodonTsvPath(actionabilityDir, refGenomeVersion);
        List<KnownCodon> knownCodons = KnownCodonFile.read(knownCodonTsv);
        LOGGER.info(" Loaded {} known codon from {}", knownCodons.size(), knownCodonTsv);

        String knownExonTsv = KnownExonFile.knownExonTsvPath(actionabilityDir, refGenomeVersion);
        List<KnownExon> knownExons = KnownExonFile.read(knownExonTsv);
        LOGGER.info(" Loaded {} known exons from {}", knownExons.size(), knownExonTsv);

        String knownCopyNumberTsv = KnownCopyNumberFile.knownCopyNumberTsvPath(actionabilityDir, refGenomeVersion);
        List<KnownCopyNumber> knownCopyNumbers = KnownCopyNumberFile.read(knownCopyNumberTsv);
        LOGGER.info(" Loaded {} known gene copy numbers from {}", knownCopyNumbers.size(), knownCopyNumberTsv);

        String knownFusionTsv = KnownFusionPairFile.knownFusionPairTsvPath(actionabilityDir, refGenomeVersion);
        List<KnownFusionPair> knownFusionPairs = KnownFusionPairFile.read(knownFusionTsv);
        LOGGER.info(" Loaded {} known fusions from {}", knownFusionPairs.size(), knownFusionTsv);

        return ImmutableKnownEvents.builder()
                .knownHotspots(hotspots)
                .knownCodons(knownCodons)
                .knownExons(knownExons)
                .knownCopyNumbers(knownCopyNumbers)
                .knownFusionPairs(knownFusionPairs)
                .build();
    }
}
