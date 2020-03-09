package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class TransvarTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(TransvarTestApplication.class);

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        String refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

        Transvar transvar = Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile);

        extractAndPrintHotspots(transvar, "MTOR", "L2230V");
        extractAndPrintHotspots(transvar, "ALK", "L1152R");
        extractAndPrintHotspots(transvar, "AKT1", "E17K");
        extractAndPrintHotspots(transvar, "ABL1", "F486S");
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @NotNull String proteinAnnotation)
            throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(gene,  proteinAnnotation);

        LOGGER.info("Printing hotspots for '{}:p.{}'", gene, proteinAnnotation);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
