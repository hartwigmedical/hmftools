package com.hartwig.hmftools.serve.transvar.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
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

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String refGenomeFastaFile37 = System.getProperty("user.home") + "/hmf/refgenomes/grch37/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        Transvar transvar37 = Transvar.withRefGenome(RefGenomeVersion.V37, refGenomeFastaFile37, HmfGenePanelSupplier.allGenesMap37());

        extractAndPrintHotspots(transvar37, "FGFR3", "ENST00000440486", "K650Q");

        String refGenomeFastaFile38 =
                System.getProperty("user.home") + "/hmf/refgenomes/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
        Transvar transvar38 = Transvar.withRefGenome(RefGenomeVersion.V38, refGenomeFastaFile38, HmfGenePanelSupplier.allGenesMap38());

        extractAndPrintHotspots(transvar38, "BRAF", "ENST00000288602", "V600E");
    }

    private static void extractAndPrintHotspots(@NotNull Transvar transvar, @NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation) {
        List<VariantHotspot> hotspots = transvar.resolve(gene, specificTranscript, proteinAnnotation);

        LOGGER.info("Printing hotspots for '{}:p.{}' on transcript {}", gene, proteinAnnotation, specificTranscript);
        for (VariantHotspot hotspot : hotspots) {
            LOGGER.info(" {}", hotspot);
        }
    }
}
