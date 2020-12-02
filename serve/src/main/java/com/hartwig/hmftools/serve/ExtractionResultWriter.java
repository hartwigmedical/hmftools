package com.hartwig.hmftools.serve;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;
import com.hartwig.hmftools.serve.hotspot.HotspotVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ExtractionResultWriter {

    private static final Logger LOGGER = LogManager.getLogger(ExtractionResultWriter.class);

    private static final String KNOWN_HOTSPOT_VCF = "KnownHotspots.SERVE.vcf";
    private static final String KNOWN_COPY_NUMBER_TSV = "KnownCopyNumbers.SERVE.tsv";
    private static final String KNOWN_FUSION_PAIR_TSV = "KnownFusionPairs.SERVE.tsv";

    private ExtractionResultWriter() {
    }

    public static void write(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion, @NotNull ExtractionResult result)
            throws IOException {
        LOGGER.info("Writing SERVE output to {}", outputDir);

        String hotspotVcfPath = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_HOTSPOT_VCF);
        LOGGER.info(" Writing {} known hotspots to {}", result.knownHotspots().size(), hotspotVcfPath);
        HotspotVCF.write(hotspotVcfPath, result.knownHotspots());

        String copyNumberTsv = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_COPY_NUMBER_TSV);
        LOGGER.info(" Writing {} known copy numbers to {}", result.knownCopyNumbers().size(), copyNumberTsv);
        // TODO Write Copy Numbers

        String fusionPairTsv = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_FUSION_PAIR_TSV);
        LOGGER.info(" Writing {} known fusion pairs to {}", result.knownFusionPairs().size(), fusionPairTsv);
        // TODO Write Fusions

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

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable signatures to {}", result.actionableSignatures().size(), actionableSignatureTsv);
        ActionableSignatureFile.write(actionableSignatureTsv, result.actionableSignatures());
    }
}
