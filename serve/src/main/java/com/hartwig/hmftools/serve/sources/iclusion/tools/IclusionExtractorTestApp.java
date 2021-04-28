package com.hartwig.hmftools.serve.sources.iclusion.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.ServeLocalConfigProvider;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.refgenome.ImmutableRefGenomeResource;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractorFactory;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionReader;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionUtil;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class IclusionExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.INFO);

        ServeConfig config = ServeLocalConfigProvider.create();

        Path outputPath = new File(config.outputDir()).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.info("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        DoidLookup doidLookup = DoidLookupFactory.buildFromMappingTsv(config.missingDoidsMappingTsv());

        IclusionExtractor extractor = IclusionExtractorFactory.buildIclusionExtractor(IclusionClassificationConfig.build(),
                buildRefGenomeResource(config),
                doidLookup);

        List<IclusionTrial> trials = IclusionReader.readAndCurate(config.iClusionTrialTsv());
        ExtractionResult result = extractor.extract(trials);

        IclusionUtil.printIclusionResult(result);
        String iclusionMutationTsv = config.outputDir() + File.separator + "IclusionMutations.tsv";
        IclusionUtil.writeIclusionMutationTypes(iclusionMutationTsv, trials);

        new ExtractionResultWriter(config.outputDir(), RefGenomeVersion.V37).write(result);
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource(@NotNull ServeConfig config) throws IOException {
        LOGGER.info("Reading driver genes from {}", config.driverGene37Tsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGene37Tsv());
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading known fusions from {}", config.knownFusion37File());
        KnownFusionCache fusionCache = new KnownFusionCache();
        if (!fusionCache.loadFile(config.knownFusion37File())) {
            throw new IllegalStateException("Could not load known fusion cache from " + config.knownFusion37File());
        }
        LOGGER.info(" Read {} known fusions", fusionCache.getData().size());

        return ImmutableRefGenomeResource.builder()
                .refSequence(new IndexedFastaSequenceFile(new File(config.refGenome37FastaFile())))
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap37())
                .proteinResolver(ProteinResolverFactory.dummy())
                .build();
    }
}
