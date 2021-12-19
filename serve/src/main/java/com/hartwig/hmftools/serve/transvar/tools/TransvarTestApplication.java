package com.hartwig.hmftools.serve.transvar.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.Transvar;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TransvarTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(TransvarTestApplication.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        Options options = TransvarTestAppConfig.createOptions();

        TransvarTestAppConfig config = null;
        try {
            config = TransvarTestAppConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("TransvarTestApp", options);
            System.exit(1);
        }

        // TODO: Load ensembl data cache correctly for 37
        Transvar transvar37 = Transvar.withRefGenome(RefGenomeVersion.V37,
                config.refGenome37FastaFile(),
                new EnsemblDataCache(Strings.EMPTY, RefGenomeVersion.V37));

        extractAndPrintHotspots(transvar37, config.gene37(), config.transcript37(), config.protein37());

        // TODO: Load ensembl data cache correctly for 38
        Transvar transvar38 = Transvar.withRefGenome(RefGenomeVersion.V38,
                config.refGenome38FastaFile(),
                new EnsemblDataCache(Strings.EMPTY, RefGenomeVersion.V38));

        extractAndPrintHotspots(transvar38, config.gene38(), config.transcript38(), config.protein38());
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
