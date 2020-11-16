package com.hartwig.hmftools.serve.sources.vicc;

import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.cli.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ViccTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccTestApplication.class);

    private static final Set<ViccSource> VICC_SOURCES_TO_INCLUDE =
            Sets.newHashSet(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);
    private static final Integer MAX_VICC_ENTRIES = null;

    public static void main(String[] args) throws IOException, ParseException {
        Configurator.setRootLevel(Level.DEBUG);

        Options options = ViccTestApplication.createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String viccJsonPath;
        String outputDir;
        ProteinResolver proteinResolver;
        List<DriverGene> driverGenes;

        if (hostname.toLowerCase().contains("datastore")) {
            viccJsonPath = "/data/common/dbs/serve/vicc/all.json";
            outputDir = System.getProperty("user.home") + "/tmp";
            proteinResolver = ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.HG19,
                    "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta");
            driverGenes = DriverGenePanelConfig.driverGenes(cmd);

        } else {
            proteinResolver = ProteinResolverFactory.dummy();
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/serve/vicc/all.json";
            outputDir = System.getProperty("user.home") + "/hmf/tmp/serve";
            driverGenes = DriverGeneFile.read(System.getProperty("user.home") + "/hmf/projects/driverGenePanel/DriverGenePanel.hg19.tsv");
        }

        Path outputPath = new File(outputDir).toPath();
        if (!Files.exists(outputPath)) {
            LOGGER.debug("Creating {} directory for writing SERVE output", outputPath.toString());
            Files.createDirectory(outputPath);
        }

        String viccFeatureTsv = outputDir + "/viccFeatures.tsv";
        String viccFeatureInterpretationTsv = outputDir + "/viccFeaturesInterpretation.tsv";

        LOGGER.debug("Configured '{}' as the VICC json path", viccJsonPath);
        LOGGER.debug("Configured '{}' as the VICC feature output TSV", viccFeatureTsv);

        List<ViccEntry> viccEntries = ViccReader.readAndCurateRelevantEntries(viccJsonPath, VICC_SOURCES_TO_INCLUDE, MAX_VICC_ENTRIES);
        ViccExtractor viccExtractor =
                ViccExtractorFactory.buildViccExtractorWithInterpretationTsv(proteinResolver, viccFeatureInterpretationTsv);

        ViccExtractionOutput viccExtractionOutput = viccExtractor.extractFromViccEntries(viccEntries, driverGenes);

        ViccUtil.writeFeatures(viccFeatureTsv, viccEntries);
        ViccUtil.writeActionability(outputDir, viccExtractionOutput);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);

        return options;
    }
}
