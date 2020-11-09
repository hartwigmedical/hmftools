package com.hartwig.hmftools.serve.transvar;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TransvarTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(TransvarTestApplication.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        String refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

        Transvar transvar = Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile);

        String gene = "KIT";
        String transcript = null;
        String proteinImpact = "T574_D575insTQLPYD";

        extractAndPrintHotspots(transvar, gene, transcript, proteinImpact);
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation) {
        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(gene, specificTranscript, proteinAnnotation);

        LOGGER.info("Printing hotspots for '{}:p.{}' on transcript {}", gene, proteinAnnotation, specificTranscript);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
