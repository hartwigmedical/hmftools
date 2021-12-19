package com.hartwig.hmftools.serve.refgenome;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class RefGenomeManagerFactory {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeManagerFactory.class);

    private RefGenomeManagerFactory() {
    }

    @NotNull
    public static RefGenomeManager createFromServeConfig(@NotNull ServeConfig config) throws IOException {
        Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap = Maps.newHashMap();
        refGenomeResourceMap.put(RefGenomeVersion.V37, buildRefGenomeResource37(config));
        refGenomeResourceMap.put(RefGenomeVersion.V38, buildRefGenomeResource38(config));
        return new RefGenomeManager(refGenomeResourceMap);
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource37(@NotNull ServeConfig config) throws IOException {
        String fastaFile37 = config.refGenome37FastaFile();
        LOGGER.info("Creating ref genome resource for V37 using fasta {}", fastaFile37);
        Map<String, HmfTranscriptRegion> transcriptMap37 = HmfGenePanelSupplier.allGenesMap37();
        ProteinResolver proteinResolver37 = config.skipHotspotResolving()
                ? ProteinResolverFactory.dummy()
                : ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.V37, fastaFile37, transcriptMap37);

        return ImmutableRefGenomeResource.builder()
                .refSequence(new IndexedFastaSequenceFile(new File(fastaFile37)))
                .driverGenes(readDriverGenesFromFile(config.driverGene37Tsv()))
                .knownFusionCache(buildKnownFusionCacheFromFile(config.knownFusion37File()))
                .ensemblDataCache(loadEnsemblDataCache(RefGenomeVersion.V37, config.ensemblDataDir37()))
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V38, config.refGenome37To38Chain())
                .proteinResolver(proteinResolver37)
                .build();
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource38(@NotNull ServeConfig config) throws IOException {
        String fastaFile38 = config.refGenome38FastaFile();
        LOGGER.info("Creating ref genome resource for V38 using fasta {}", fastaFile38);
        Map<String, HmfTranscriptRegion> transcriptMap38 = HmfGenePanelSupplier.allGenesMap38();
        ProteinResolver proteinResolver38 = config.skipHotspotResolving()
                ? ProteinResolverFactory.dummy()
                : ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.V38, fastaFile38, transcriptMap38);

        return ImmutableRefGenomeResource.builder()
                .refSequence(new IndexedFastaSequenceFile(new File(fastaFile38)))
                .driverGenes(readDriverGenesFromFile(config.driverGene38Tsv()))
                .knownFusionCache(buildKnownFusionCacheFromFile(config.knownFusion38File()))
                .ensemblDataCache(loadEnsemblDataCache(RefGenomeVersion.V38, config.ensemblDataDir38()))
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V37, config.refGenome38To37Chain())
                .proteinResolver(proteinResolver38)
                .build();
    }

    @NotNull
    private static List<DriverGene> readDriverGenesFromFile(@NotNull String driverGeneTsv) throws IOException {
        LOGGER.info(" Reading driver genes from {}", driverGeneTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsv);
        LOGGER.info("  Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    private static KnownFusionCache buildKnownFusionCacheFromFile(@NotNull String knownFusionFile) throws IOException {
        LOGGER.info(" Reading known fusions from {}", knownFusionFile);
        KnownFusionCache cache = new KnownFusionCache();
        if (!cache.loadFile(knownFusionFile)) {
            throw new IOException("Could not load known fusions from " + knownFusionFile);
        }
        LOGGER.info("  Read {} known fusion entries", cache.getData().size());
        return cache;
    }

    @NotNull
    private static EnsemblDataCache loadEnsemblDataCache(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String ensemblDataDir)
            throws IOException {
        LOGGER.info(" Reading ensembl data cache from {}", ensemblDataDir);
        EnsemblDataCache ensemblDataCache = EnsemblDataCacheLoader.load(ensemblDataDir, refGenomeVersion);
        int geneCount = 0;
        for (List<GeneData> genesPerChromosome : ensemblDataCache.getChrGeneDataMap().values()) {
            geneCount += genesPerChromosome.size();
        }
        LOGGER.info("  Loaded entries for {} genes", geneCount);
        return ensemblDataCache;
    }
}
