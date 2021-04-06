package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.IOException;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        String ckbDir = "/data/common/dbs/ckb/210326_flex_dump";
        String outputDir = System.getProperty("user.home") + "/tmp/serve_ckb";
        String eventsTsv = outputDir + "/CkbEvents.tsv";

        // Read and curate CKB datamodel
        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        List<CkbEntry> allCkbEntries = JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
        List<CkbEntry> filteredAndcurateCkbEntries = CkbReader.filterAndCurateRelevantEntries(allCkbEntries);

        // Read required data
        String missingDoidMappingTsv = "/data/common/dbs/serve/curation/missing_doids_mapping.tsv";
        DoidLookup doidLookup = DoidLookupFactory.buildFromConfigTsv(missingDoidMappingTsv);

        String driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.hg38.tsv";
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        KnownFusionCache fusionCache = new KnownFusionCache();
        String knownFusionFilePath = "/data/common/dbs/fusions/known_fusion_data.38_v3.csv";
        if (!fusionCache.loadFile(knownFusionFilePath)) {
            throw new IllegalStateException("Could not load known fusion cache from " + knownFusionFilePath);
        }
        LOGGER.debug(" Read {} known fusions from {}", fusionCache.getData().size(), knownFusionFilePath);

        String fastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
        RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38;
        Map<String, HmfTranscriptRegion> allGenesMap = HmfGenePanelSupplier.allGenesMap38();
        ProteinResolver proteinResolver = ProteinResolverFactory.transvarWithRefGenome(refGenomeVersion, fastaFile, allGenesMap);
        RefGenomeResource refGenomeResource = ImmutableRefGenomeResource.builder()
                .fastaFile(fastaFile)
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .canonicalTranscriptPerGeneMap(allGenesMap)
                .proteinResolver(proteinResolver)
                .build();

        EventClassifierConfig config = CkbClassificationConfig.build();
        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(config, refGenomeResource, doidLookup);

        ExtractionResult result = extractor.extract(filteredAndcurateCkbEntries);

        CkbUtils.writeEventsToTsv(eventsTsv, filteredAndcurateCkbEntries);
        CkbUtils.printExtractionResults(result);

        new ExtractionResultWriter(outputDir, refGenomeVersion).write(result);
    }
}