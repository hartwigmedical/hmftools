package com.hartwig.hmftools.serve.sources.vicc.tools;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;
import com.hartwig.hmftools.serve.sources.vicc.ViccUtil;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class ViccExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractorTestApp.class);

    private static final Set<ViccSource> VICC_SOURCES_TO_INCLUDE =
            Sets.newHashSet(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);
    private static final Integer MAX_VICC_ENTRIES = null;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String viccJsonPath;
        String driverGeneTsvPath;
        String outputDir;
        ProteinResolver proteinResolver;

        if (hostname.toLowerCase().contains("datastore")) {
            viccJsonPath = "/data/common/dbs/serve/vicc/all.json";
            driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.hg19.tsv";
            outputDir = System.getProperty("user.home") + "/tmp";
            proteinResolver = ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.HG19,
                    "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta");
        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/serve/static_sources/vicc/all.json";
            driverGeneTsvPath = System.getProperty("user.home") + "/hmf/projects/driverGenePanel/DriverGenePanel.hg19.tsv";
            outputDir = System.getProperty("user.home") + "/hmf/tmp/serve";
            proteinResolver = ProteinResolverFactory.dummy();
        }

        Path outputPath = new File(outputDir).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.debug("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        String featureTsv = outputDir + "/viccFeatures.tsv";
        String featureInterpretationTsv = outputDir + "/viccFeatureInterpretation.tsv";

        LOGGER.debug("Configured '{}' as the VICC json path", viccJsonPath);
        LOGGER.debug("Configured '{}' as the VICC feature output TSV path", featureTsv);
        LOGGER.debug("Configured '{}' as the VICC feature interpretation output TSV path", featureInterpretationTsv);
        LOGGER.debug("Configured '{}' as the driver gene TSV path", driverGeneTsvPath);

        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        Map<String, HmfTranscriptRegion> allGenesMap = HmfGenePanelSupplier.allGenesMap37();
        List<ViccEntry> entries = ViccReader.readAndCurateRelevantEntries(viccJsonPath, VICC_SOURCES_TO_INCLUDE, MAX_VICC_ENTRIES);
        ViccExtractor viccExtractor = ViccExtractorFactory.buildViccExtractorWithInterpretationTsv(proteinResolver,
                driverGenes,
                allGenesMap,
                featureInterpretationTsv);

        ExtractionOutput extractionOutput = viccExtractor.extractFromViccEntries(entries);

        ViccUtil.writeFeatures(featureTsv, entries);
        ViccUtil.writeActionability(outputDir, extractionOutput);
    }
}
