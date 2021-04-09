package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.ckb.JsonDatabaseToCkbEntryConverter;
import com.hartwig.hmftools.ckb.classification.CkbClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.CkbJsonReader;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.refseq.RefSeq;
import com.hartwig.hmftools.common.refseq.RefSeqFile;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.refgenome.ImmutableRefGenomeResource;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractor;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractorFactory;
import com.hartwig.hmftools.serve.sources.ckb.CkbReader;
import com.hartwig.hmftools.serve.sources.ckb.CkbUtils;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);
        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String ckbDir;
        String outputDir;
        String missingDoidMappingTsv;
        String driverGeneTsvPath;
        String knownFusionFilePath;
        String fastaFile;
        ProteinResolver proteinResolver;
        String eventsTsv;
        String refSeqMatch;

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38;
        Map<String, HmfTranscriptRegion> allGenesMap = HmfGenePanelSupplier.allGenesMap38();

        if (hostname.toLowerCase().contains("datastore")) {
            ckbDir = "/data/common/dbs/ckb/210402_flex_dump";
            outputDir = System.getProperty("user.home") + "/tmp/serve_ckb";
            eventsTsv = outputDir + "/CkbEvents.tsv";
            missingDoidMappingTsv = "/data/common/dbs/serve/curation/missing_doids_mapping.tsv";
            driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.38.tsv";
            knownFusionFilePath = "/data/common/dbs/fusions/known_fusion_data.38_v3.csv";
            fastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
            proteinResolver = ProteinResolverFactory.dummy();
            refSeqMatch = "/data/common/dbs/serve/static_sources/refseq/refseq_to_canonicalTranscript.tsv";
        } else {
            ckbDir = System.getProperty("user.home") + "/hmf/projects/serve/ckb";
            outputDir = System.getProperty("user.home") + "/tmp/serve_ckb";
            eventsTsv = outputDir + "/CkbEvents.tsv";
            missingDoidMappingTsv = System.getProperty("user.home") + "/hmf/projects/serve/curation/missing_doids_mapping.tsv";
            driverGeneTsvPath = System.getProperty("user.home") + "/hmf/projects/driverGenePanel/DriverGenePanel.38.tsv";
            knownFusionFilePath = System.getProperty("user.home") + "/hmf/projects/fusions/known_fusion_data.38_v3.csv";
            fastaFile = Strings.EMPTY;
            proteinResolver = ProteinResolverFactory.dummy();
            refSeqMatch = System.getProperty("user.home") + "/hmf/projects/serve/static_sources/refseq/refseq_to_canonicalTranscript.tsv";
        }

        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        List<CkbEntry> allCkbEntries = JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
        List<CkbEntry> filteredAndcurateCkbEntries = CkbReader.filterAndCurateRelevantEntries(allCkbEntries, 1000);

        DoidLookup doidLookup = DoidLookupFactory.buildFromConfigTsv(missingDoidMappingTsv);

        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        KnownFusionCache fusionCache = new KnownFusionCache();
        if (!fusionCache.loadFile(knownFusionFilePath)) {
            throw new IllegalStateException("Could not load known fusion cache from " + knownFusionFilePath);
        }
        LOGGER.debug(" Read {} known fusions from {}", fusionCache.getData().size(), knownFusionFilePath);

        RefGenomeResource refGenomeResource = ImmutableRefGenomeResource.builder()
                .fastaFile(fastaFile)
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .canonicalTranscriptPerGeneMap(allGenesMap)
                .proteinResolver(proteinResolver)
                .build();

        EventClassifierConfig config = CkbClassificationConfig.build();
        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(config, refGenomeResource, doidLookup);

        LOGGER.info("Reading ref seq matching to transcript");
        List<RefSeq> refSeqMatchFile = RefSeqFile.readingRefSeq(refSeqMatch);

        ExtractionResult result = extractor.extract(filteredAndcurateCkbEntries, refSeqMatchFile);

        CkbUtils.writeEventsToTsv(eventsTsv, filteredAndcurateCkbEntries);
        CkbUtils.printExtractionResults(result);

        new ExtractionResultWriter(outputDir, refGenomeVersion).write(result);
    }
}