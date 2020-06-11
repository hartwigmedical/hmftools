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

        // Below cases are different genomic mutations leading to identical protein impacts due to redundancy in trinucleotide encoding
        // (different bases encoding for the same AA)
        // Eg, While the mutations that delete codons 842-845 from PDGFRA are genomically different from 843-846, in terms of AA they will
        // have the same impact due to redundancy in trinucleotide encoding.
//        extractAndPrintHotspots(transvar, "PDGFRA", "ENST00000257290", "D842_H845del");
//        extractAndPrintHotspots(transvar, "PDGFRA", "ENST00000257290", "I843_D846del");
//
//        extractAndPrintHotspots(transvar, "KIT", "ENST00000288135", "K550_W557del");
//        extractAndPrintHotspots(transvar, "KIT", "ENST00000288135", "P551_K558del");
//
//        extractAndPrintHotspots(transvar, "KIT", null, "V555_V559del");
//        extractAndPrintHotspots(transvar, "KIT", null, "Q556_V560del");

//        extractAndPrintHotspots(transvar, "KIT", null, "P577_W582delinsPYD");
//        extractAndPrintHotspots(transvar, "KIT", null, "H580_W582del");

        // Below variant is annotated as p.DI842IM by SnpEff which is functionally identical to below.
//        extractAndPrintHotspots(transvar, "PDGFRA", "ENST00000257290", "D842_I843delinsIM");

        extractAndPrintHotspots(transvar, "CREBBP", null, "S1680del");
        extractAndPrintHotspots(transvar, "CREBBP", null, "S1679del");




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
