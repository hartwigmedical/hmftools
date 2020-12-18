package com.hartwig.hmftools.serve;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.checkertool.CheckCodonRanges;
import com.hartwig.hmftools.serve.checkertool.CheckExons;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

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

        Map<String, HmfTranscriptRegion> allGenesMap = selectAllGeneMap(config.refGenomeVersion());

        ServeAlgo algo = new ServeAlgo(readDriverGenesFromFile(config.driverGeneTsv()),
                buildKnownFusionCacheFromFile(config.knownFusionFile()),
                allGenesMap,
                buildProteinResolver(config, allGenesMap),
                buildDoidLookup(config.missingDoidsMappingTsv()));

        ExtractionResult result = algo.run(config);

        CheckExons checkExons = new CheckExons();
        checkExons.run(result);

        CheckCodonRanges checkCodonRanges = new CheckCodonRanges();
        checkCodonRanges.run(result);

        new ExtractionResultWriter(config.outputDir(), config.refGenomeVersion()).write(result);

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
    private static Map<String, HmfTranscriptRegion> selectAllGeneMap(@NotNull RefGenomeVersion refGenomeVersion) {
        if (refGenomeVersion == RefGenomeVersion.V37) {
            return HmfGenePanelSupplier.allGenesMap37();
        }

        LOGGER.warn("Ref genome version not supported: '{}'. Reverting to V37 Gene Map.", refGenomeVersion);
        return HmfGenePanelSupplier.allGenesMap37();
    }

    @NotNull
    private static ProteinResolver buildProteinResolver(@NotNull ServeConfig config,
            @NotNull Map<String, HmfTranscriptRegion> transcriptsPerGeneMap) throws FileNotFoundException {
        if (config.skipHotspotResolving()) {
            LOGGER.info("Creating dummy protein resolver");
            return ProteinResolverFactory.dummy();
        } else {
            LOGGER.info("Creating transvar protein resolver with ref genome version '{}' using FASTA {}",
                    config.refGenomeVersion(),
                    config.refGenomeFastaFile());
            return ProteinResolverFactory.transvarWithRefGenome(config.refGenomeVersion(),
                    config.refGenomeFastaFile(),
                    transcriptsPerGeneMap);
        }
    }

    @NotNull
    private static DoidLookup buildDoidLookup(@NotNull String missingDoidsMappingTsv) throws IOException {
        LOGGER.info("Creating missing doid lookup mapping from {}", missingDoidsMappingTsv);
        return DoidLookupFactory.buildFromConfigTsv(missingDoidsMappingTsv);
    }
}
