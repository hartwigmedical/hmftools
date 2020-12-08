package com.hartwig.hmftools.serve.extraction;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumberFile;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPairFile;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspotVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ExtractionResultWriter {

    private static final Logger LOGGER = LogManager.getLogger(ExtractionResultWriter.class);

    private static final String KNOWN_HOTSPOT_VCF = "KnownHotspots.SERVE.vcf.gz";
    private static final String KNOWN_COPY_NUMBER_TSV = "KnownCopyNumbers.SERVE.tsv";
    private static final String KNOWN_FUSION_PAIR_TSV = "KnownFusionPairs.SERVE.tsv";

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

        String hotspotVcfPath = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_HOTSPOT_VCF);
        LOGGER.info(" Writing {} known hotspots to {}", result.knownHotspots().size(), hotspotVcfPath);
        KnownHotspotVCF.write(hotspotVcfPath, result.knownHotspots());

        String copyNumberTsv = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_COPY_NUMBER_TSV);
        LOGGER.info(" Writing {} known copy numbers to {}", result.knownCopyNumbers().size(), copyNumberTsv);
        KnownCopyNumberFile.write(copyNumberTsv, result.knownCopyNumbers());

        String fusionPairTsv = refGenomeVersion.makeVersioned(outputDir + File.separator + KNOWN_FUSION_PAIR_TSV);
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

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(outputDir, refGenomeVersion);
        LOGGER.info(" Writing {} actionable signatures to {}", result.actionableSignatures().size(), actionableSignatureTsv);
        ActionableSignatureFile.write(actionableSignatureTsv, result.actionableSignatures());
    }
}
