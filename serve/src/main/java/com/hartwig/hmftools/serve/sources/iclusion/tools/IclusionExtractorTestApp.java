package com.hartwig.hmftools.serve.sources.iclusion.tools;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractorFactory;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionReader;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionUtil;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class IclusionExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String iclusionTrialTsv;
        String driverGeneTsvPath;
        String missingDoidMappingTsv;
        String outputDir;
        ProteinResolver proteinResolver;

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.V37;
        Map<String, HmfTranscriptRegion> allGenesMap = HmfGenePanelSupplier.allGenesMap37();
        if (hostname.toLowerCase().contains("datastore")) {
            iclusionTrialTsv = "/data/common/dbs/iclusion/iclusion_trials_prod.tsv";
            driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.hg19.tsv";
            missingDoidMappingTsv = "/data/common/dbs/serve/curation/missing_doids_mapping.tsv";
            outputDir = System.getProperty("user.home") + "/tmp";
            proteinResolver = ProteinResolverFactory.transvarWithRefGenome(refGenomeVersion,
                    "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta",
                    allGenesMap);
        } else {
            iclusionTrialTsv = System.getProperty("user.home") + "/hmf/projects/serve/iclusion/iclusion_trials_prod.tsv";
            driverGeneTsvPath = System.getProperty("user.home") + "/hmf/projects/driverGenePanel/DriverGenePanel.hg19.tsv";
            missingDoidMappingTsv = System.getProperty("user.home") + "/hmf/projects/serve/curation/missing_doids_mapping.tsv";
            outputDir = System.getProperty("user.home") + "/hmf/tmp/serve";
            proteinResolver = ProteinResolverFactory.dummy();
        }

        Path outputPath = new File(outputDir).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.debug("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        String iclusionMutationTsv = outputDir + "/iclusionMutations.tsv";

        LOGGER.debug("Configured '{}' as the iClusion trial TSV path", iclusionTrialTsv);
        LOGGER.debug("Configured '{}' as the driver gene TSV path", driverGeneTsvPath);
        LOGGER.debug("Configured '{}' as the missing DOID mapping TSV path", missingDoidMappingTsv);
        LOGGER.debug("Configured '{}' as the iClusion mutation TSV path", iclusionMutationTsv);

        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        List<IclusionTrial> trials = IclusionReader.readAndCurate(iclusionTrialTsv);

        DoidLookup doidLookup = DoidLookupFactory.buildFromConfigTsv(missingDoidMappingTsv);

        EventClassifierConfig config = IclusionClassificationConfig.build();
        IclusionExtractor extractor =
                IclusionExtractorFactory.buildIclusionExtractor(config, proteinResolver, driverGenes, allGenesMap, doidLookup);
        ExtractionResult result = extractor.extract(trials);

        IclusionUtil.printIclusionResult(result);
        IclusionUtil.writeIclusionMutationTypes(iclusionMutationTsv, trials);

        new ExtractionResultWriter(outputDir, refGenomeVersion).write(result);
    }

}
