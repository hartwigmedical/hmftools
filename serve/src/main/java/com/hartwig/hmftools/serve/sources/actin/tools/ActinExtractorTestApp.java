package com.hartwig.hmftools.serve.sources.actin.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.refgenome.ImmutableRefGenomeResource;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.actin.ActinExtractor;
import com.hartwig.hmftools.serve.sources.actin.ActinExtractorFactory;
import com.hartwig.hmftools.serve.sources.actin.ActinReader;
import com.hartwig.hmftools.serve.sources.actin.ActinUtils;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinClassificationConfig;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ActinExtractorTestApp {
    private static final Logger LOGGER = LogManager.getLogger(ActinExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        ServeConfig config = ServeLocalConfigProvider.create();

        Path outputPath = new File(config.outputDir()).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.info("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        RefGenomeResource refGenomeResource = buildRefGenomeResource(config);
        ActinExtractor extractor = ActinExtractorFactory.buildActinExtractor(ActinClassificationConfig.build(), refGenomeResource);

        List<ActinEntry> actinEntries = ActinReader.read(config.actinTrialTsv());
        ExtractionResult result = extractor.extract(actinEntries);

        String eventsTsv = config.outputDir() + File.separator + "ActinEventClassification.tsv";
        ActinUtils.writeActinEventTypes(eventsTsv, actinEntries);

        new ExtractionResultWriter(config.outputDir(), Knowledgebase.ACTIN.refGenomeVersion(), refGenomeResource.refSequence()).write(result);
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource(@NotNull ServeConfig config) throws IOException {
        LOGGER.info("Reading driver genes from {}", config.driverGene38Tsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGene38Tsv());
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading known fusions from {}", config.knownFusion38File());
        KnownFusionCache fusionCache = new KnownFusionCache();
        if (!fusionCache.loadFile(config.knownFusion38File())) {
            throw new IllegalStateException("Could not load known fusion cache from " + config.knownFusion38File());
        }
        LOGGER.info(" Read {} known fusions", fusionCache.getData().size());

        return ImmutableRefGenomeResource.builder()
                .refSequence(new IndexedFastaSequenceFile(new File(config.refGenome38FastaFile())))
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap38())
                .proteinResolver(ProteinResolverFactory.dummy())
                .build();
    }


}
