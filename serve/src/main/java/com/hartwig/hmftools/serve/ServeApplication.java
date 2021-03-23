package com.hartwig.hmftools.serve;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.refgenome.RefGenomeManager;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ServeApplication {

    private static final Logger LOGGER = LogManager.getLogger(ServeApplication.class);

    private static final String VERSION = ServeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE v{}", VERSION);

        Options options = ServeConfig.createOptions();

        ServeConfig config = null;
        try {
            config = ServeConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("SERVE", options);
            System.exit(1);
        }

        RefGenomeManager refGenomeManager = RefGenomeManager.buildFromServeConfig(config);

        ServeAlgo algo = new ServeAlgo(refGenomeManager,
                readDriverGenesFromFile(config.driverGeneTsv()),
                buildKnownFusionCacheFromFile(config.knownFusionFile()),
                buildDoidLookup(config.missingDoidsMappingTsv()));

        ExtractionResult result = algo.run(config);

        new ExtractionResultWriter(config.outputDir(), RefGenomeVersion.V37).write(result);

        LOGGER.info("Done!");
    }

    @NotNull
    private static List<DriverGene> readDriverGenesFromFile(@NotNull String driverGeneTsv) throws IOException {
        LOGGER.info("Reading driver genes from {}", driverGeneTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsv);
        LOGGER.info(" Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    private static KnownFusionCache buildKnownFusionCacheFromFile(@NotNull String knownFusionFile) throws IOException {
        LOGGER.info("Reading known fusions from {}", knownFusionFile);
        KnownFusionCache cache = new KnownFusionCache();
        if (!cache.loadFile(knownFusionFile)) {
            throw new IOException("Could not load known fusions from " + knownFusionFile);
        }
        LOGGER.info(" Read {} known fusion entries", cache.getData().size());
        return cache;
    }

    @NotNull
    private static DoidLookup buildDoidLookup(@NotNull String missingDoidsMappingTsv) throws IOException {
        LOGGER.info("Creating missing doid lookup mapping from {}", missingDoidsMappingTsv);
        return DoidLookupFactory.buildFromConfigTsv(missingDoidsMappingTsv);
    }
}
