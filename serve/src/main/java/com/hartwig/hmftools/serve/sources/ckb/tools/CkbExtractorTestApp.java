package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.ckb.classification.CkbClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.refgenome.EnsemblDataCacheLoader;
import com.hartwig.hmftools.serve.refgenome.ImmutableRefGenomeResource;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractor;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractorFactory;
import com.hartwig.hmftools.serve.sources.ckb.CkbReader;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        ServeConfig config = ServeLocalConfigProvider.create();

        Path outputPath = new File(config.outputDir()).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.info("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        RefGenomeResource refGenomeResource = buildRefGenomeResource(config);
        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(CkbClassificationConfig.build(), refGenomeResource);

        List<CkbEntry> entries = CkbReader.readAndCurate(config.ckbDir(), config.ckbFilterTsv());
        ExtractionResult result = extractor.extract(entries, Maps.newHashMap());

        String eventsTsv = config.outputDir() + File.separator + "CkbEventClassification.tsv";
        CkbUtil.writeEventsToTsv(eventsTsv, entries);
        CkbUtil.printExtractionResults(result);

        new ExtractionResultWriter(config.outputDir(), Knowledgebase.CKB.refGenomeVersion(), refGenomeResource.refSequence()).write(result);
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource(@NotNull ServeConfig config) throws IOException {
        LOGGER.info("Reading driver genes from {}", config.driverGene38Tsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGene38Tsv());
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading known fusions from {}", config.knownFusion38File());
        KnownFusionCache fusionCache = new KnownFusionCache();
        if (!fusionCache.loadFile(config.knownFusion38File())) {
            throw new IOException("Could not load known fusion cache from " + config.knownFusion38File());
        }
        LOGGER.info(" Read {} known fusions", fusionCache.getData().size());

        LOGGER.info(" Reading ensembl data cache from {}", config.ensemblDataDir38());
        EnsemblDataCache ensemblDataCache = EnsemblDataCacheLoader.load(config.ensemblDataDir38(), RefGenomeVersion.V38);
        LOGGER.info("  Loaded ensembl data cache from {}", ensemblDataCache);

        return ImmutableRefGenomeResource.builder()
                .refSequence(new IndexedFastaSequenceFile(new File(config.refGenome38FastaFile())))
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .ensemblDataCache(ensemblDataCache)
                .proteinResolver(ProteinResolverFactory.dummy())
                .build();
    }
}