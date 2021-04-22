package com.hartwig.hmftools.serve.extraction;

import java.io.IOException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristicFile;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableSignatureFile;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodonFile;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumberFile;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPairFile;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspotFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ExtractionResultWriter {

    private static final Logger LOGGER = LogManager.getLogger(ExtractionResultWriter.class);

    @NotNull
    private final String outputDir;
    @NotNull
    private final RefGenomeVersion refGenomeVersion;

    public ExtractionResultWriter(@NotNull final String outputDir, @NotNull final RefGenomeVersion refGenomeVersion) {
        this.outputDir = outputDir;
        this.refGenomeVersion = refGenomeVersion;
    }

    public void write(@NotNull ExtractionResult result) throws IOException {
        LOGGER.info("Writing SERVE output to {}", outputDir);

        String hotspotVcf = KnownHotspotFile.knownHotspotVcfPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} known hotspots to {}", result.knownHotspots().size(), hotspotVcf);
        KnownHotspotFile.write(hotspotVcf, result.knownHotspots());

        String codonTsv = KnownCodonFile.knownCodonTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} known codons to {}", result.knownCodons().size(), codonTsv);
        KnownCodonFile.write(codonTsv, result.knownCodons());

        String exonTsv = KnownExonFile.knownExonTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} known exons to {}", result.knownExons().size(), exonTsv);
        KnownExonFile.write(exonTsv, result.knownExons());

        String copyNumberTsv = KnownCopyNumberFile.knownCopyNumberTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} known copy numbers to {}", result.knownCopyNumbers().size(), copyNumberTsv);
        KnownCopyNumberFile.write(copyNumberTsv, result.knownCopyNumbers());

        String fusionPairTsv = KnownFusionPairFile.knownFusionPairTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} known fusion pairs to {}", result.knownFusionPairs().size(), fusionPairTsv);
        KnownFusionPairFile.write(fusionPairTsv, result.knownFusionPairs());

        String actionableHotspotTsv = ActionableHotspotFile.actionableHotspotTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable hotspots to {}", result.actionableHotspots().size(), actionableHotspotTsv);
        ActionableHotspotFile.write(actionableHotspotTsv, result.actionableHotspots());

        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable ranges to {}", result.actionableRanges().size(), actionableRangeTsv);
        ActionableRangeFile.write(actionableRangeTsv, result.actionableRanges());

        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable genes to {}", result.actionableGenes().size(), actionableGeneTsv);
        ActionableGeneFile.write(actionableGeneTsv, result.actionableGenes());

        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable fusions to {}", result.actionableFusions().size(), actionableFusionTsv);
        ActionableFusionFile.write(actionableFusionTsv, result.actionableFusions());

        // TODO: Can be removed once PROTECT reads actionable characteristics
        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable signatures to {}", result.actionableCharacteristics().size(), actionableSignatureTsv);
        ActionableSignatureFile.write(actionableSignatureTsv, result.actionableCharacteristics());

        String actionableCharacteristicTsv = ActionableCharacteristicFile.actionableCharacteristicTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable tumor characteristics to {}",
                result.actionableCharacteristics().size(),
                actionableCharacteristicTsv);
        ActionableCharacteristicFile.write(actionableCharacteristicTsv, result.actionableCharacteristics());
    }
}
