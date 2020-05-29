package com.hartwig.hmftools.serve;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.Transvar;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TransvarTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(TransvarTestApplication.class);

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        String refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

        Transvar transvar = Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile);

        // These 2 variants are identical (the trinucleotide is repeated)
//        extractAndPrintHotspots(transvar, "KIT", null, "V560del");
//        extractAndPrintHotspots(transvar, "KIT", null, "V559del");

        // These variants lead to the same gDNA because D1739Y is not defined on the canonical transcript of BRCA1
//        extractAndPrintHotspots(transvar, "BRCA1", "ENST00000357654", "D1739Y");
//        extractAndPrintHotspots(transvar, "BRCA1", "ENST00000357654", "D1692Y");

        extractAndPrintHotspots(transvar, "EGFR", null, "L747_A750delinsP");
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @Nullable String transcriptId,
            @NotNull String proteinAnnotation) throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(gene, transcriptId, proteinAnnotation);

        LOGGER.info("Printing hotspots for '{}:p.{}'", gene, proteinAnnotation);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
