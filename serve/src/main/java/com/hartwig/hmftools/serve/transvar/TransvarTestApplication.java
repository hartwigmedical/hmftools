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

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        String refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

        Transvar transvar = Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile);


        // Leads to a warning which we can ignore.
        extractAndPrintHotspots(transvar, "BRAF", null, "T599insTT");

        // Repeat issue - These 2 variants are identical (the trinucleotide is repeated)
//        extractAndPrintHotspots(transvar, "KIT", null, "V560del");
//        extractAndPrintHotspots(transvar, "KIT", null, "V559del");
//
//        // Transcript issue - These variants lead to the same gDNA because D1739Y is not defined on the canonical transcript of BRCA1
//        extractAndPrintHotspots(transvar, "BRCA1", "ENST00000357654", "D1739Y");
//        extractAndPrintHotspots(transvar, "BRCA1", "ENST00000357654", "D1692Y");
//
//        // This variant doesn't work properly yet (complex del-then-insertion)
//        extractAndPrintHotspots(transvar, "EGFR", null, "L747_A750delinsP");
//
//        // Transcript issue
//        extractAndPrintHotspots(transvar, "FBXW7", null, "R658Q");
//        extractAndPrintHotspots(transvar, "FBXW7", null, "R482Q");
//
//        // Transcript issue
//        extractAndPrintHotspots(transvar, "FGFR2", null, "M537I");
//        extractAndPrintHotspots(transvar, "FGFR2", null, "M535I");
//
//        // Transcript issue
//        extractAndPrintHotspots(transvar, "GNAS", null, "R201H");
//        extractAndPrintHotspots(transvar, "GNAS", null, "R844H");
//
//        // Repeat issue
//        extractAndPrintHotspots(transvar, "PTEN", null, "I33del");
//        extractAndPrintHotspots(transvar, "PTEN", null, "I32del");
//
//        // Transcript issue - Complex?
//        extractAndPrintHotspots(transvar, "RUNX1", null, "R174*");
//        extractAndPrintHotspots(transvar, "RUNX1", null, "R177*");
//        extractAndPrintHotspots(transvar, "RUNX1", null, "R177Q");
//        extractAndPrintHotspots(transvar, "RUNX1", null, "R201Q");
//        extractAndPrintHotspots(transvar, "RUNX1", null, "R174Q");
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation) throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(gene, specificTranscript, proteinAnnotation);

        LOGGER.info("Printing hotspots for '{}:p.{}' on transcript {}", gene, proteinAnnotation, specificTranscript);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
