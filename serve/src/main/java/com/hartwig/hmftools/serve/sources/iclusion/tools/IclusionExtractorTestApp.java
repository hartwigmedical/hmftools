package com.hartwig.hmftools.serve.sources.iclusion.tools;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.curation.IclusionCurator;
import com.hartwig.hmftools.serve.sources.iclusion.filter.IclusionFilter;

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
        String outputDir;
        ProteinResolver proteinResolver;

        if (hostname.toLowerCase().contains("datastore")) {
            iclusionTrialTsv = "/data/common/dbs/iclusion/iclusion_trials_prod.tsv";
            driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.hg19.tsv";
            outputDir = System.getProperty("user.home") + "/tmp";
            proteinResolver = ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.HG19,
                    "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta");
        } else {
            iclusionTrialTsv = System.getProperty("user.home") + "/hmf/projects/serve/iclusion/iclusion_trials_prod.tsv";
            driverGeneTsvPath = System.getProperty("user.home") + "/hmf/projects/driverGenePanel/DriverGenePanel.hg19.tsv";
            outputDir = System.getProperty("user.home") + "/hmf/tmp/serve";
            proteinResolver = ProteinResolverFactory.dummy();
        }

        Path outputPath = new File(outputDir).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.debug("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        LOGGER.debug("Configured '{}' as the iclusion trial TSV path", iclusionTrialTsv);
        LOGGER.debug("Configured '{}' as the driver gene TSV path", driverGeneTsvPath);

        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        List<IclusionTrial> trials = IclusionTrialFile.read(iclusionTrialTsv);
        LOGGER.info(" Read {} trials from {}", trials.size(), iclusionTrialTsv);

        IclusionFilter filter = new IclusionFilter();
        List<IclusionTrial> filteredTrials = filter.run(trials);
        LOGGER.info("  {} trials remaining after filtering", filteredTrials.size());
        filter.reportUnusedFilterEntries();

        LOGGER.info(" Starting iClusion trial curation");
        IclusionCurator curator = new IclusionCurator();
        List<IclusionTrial> curatedTrials = curator.run(filteredTrials);
        curator.reportUnusedCurationEntries();

        ExtractionOutput output = new IclusionExtractor().extractFromIclusionTrials(curatedTrials);
        LOGGER.info("Generated {}", output);
    }
}
