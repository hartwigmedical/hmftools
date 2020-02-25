package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.RefVersion;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class TempRunner {

    private static final Logger LOGGER = LogManager.getLogger(TempRunner.class);

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        String refFastaPath = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        RefVersion refVersion = RefVersion.HG19;

        Transvar transvar = new Transvar(refFastaPath, refVersion);

        extractAndPrintHotspots(transvar, "MTOR", "L2230V");
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @NotNull String proteinAnnotation)
            throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(gene,  proteinAnnotation);

        LOGGER.info("Printing hotspots for {}:p.{}", gene, proteinAnnotation);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
