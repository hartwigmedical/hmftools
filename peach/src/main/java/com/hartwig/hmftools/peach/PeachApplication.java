package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.peach.data_loader.HaplotypeEventLoader;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.output.AllHaplotypeCombinationsFile;
import com.hartwig.hmftools.peach.output.BestHaplotypeCombinationsFile;
import com.hartwig.hmftools.peach.output.EventsFile;
import com.hartwig.hmftools.peach.output.EventsPerGeneFile;
import com.hartwig.hmftools.peach.output.QcStatusFile;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PeachApplication
{
    @NotNull
    private final PeachConfig config;

    public PeachApplication(@NotNull final ConfigBuilder configBuilder)
    {
        this.config = new PeachConfig(configBuilder);
    }

    public void run()
    {
        try
        {
            if(!config.isValid())
            {
                throw new IllegalArgumentException("invalid config");
            }
            PCH_LOGGER.info("running PEACH");

            PCH_LOGGER.info("creating output directory");
            createOutputDirectory();

            PCH_LOGGER.info("load haplotype config");
            HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(config.haplotypesFile);

            PCH_LOGGER.info("load events");
            Map<String, Integer> eventIdToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                    config.vcfFile, config.sampleName, haplotypePanel.getRelevantVariantPositions()
            );

            PCH_LOGGER.info("call haplotypes");
            HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
            Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = caller.getGeneToHaplotypeAnalysis(eventIdToCount);

            writeOutputFiles(eventIdToCount, geneToHaplotypeAnalysis);
        }
        catch(Exception e)
        {
            PCH_LOGGER.error("unrecoverable error encountered", e);
            System.exit(1);
        }

        PCH_LOGGER.info("finished running PEACH");
    }

    private void createOutputDirectory()
    {
        // Output directory cannot be null in valid config
        assert config.outputDir != null;
        File outputDirectory = new File(config.outputDir);
        if(!outputDirectory.exists() && !outputDirectory.mkdirs())
        {
            throw new RuntimeException(String.format("could not create output directory: %s", config.outputDir));
        }
    }

    private void writeOutputFiles(Map<String, Integer> eventIdToCount, Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis)
    {
        try
        {
            PCH_LOGGER.info("write events output file");
            EventsFile.write(config.getEventsOutputPath(), eventIdToCount);
            PCH_LOGGER.info("write events per gene output file");
            EventsPerGeneFile.write(config.getEventsPerGeneOutputPath(), geneToHaplotypeAnalysis);
            PCH_LOGGER.info("write all haplotype combinations output file");
            AllHaplotypeCombinationsFile.write(config.getAllHaplotypeCombinationsOutputPath(), geneToHaplotypeAnalysis);
            PCH_LOGGER.info("write best haplotype combination output file");
            BestHaplotypeCombinationsFile.write(config.getBestHaplotypeCombinationsOutputPath(), geneToHaplotypeAnalysis);
            PCH_LOGGER.info("write qc status output file");
            QcStatusFile.write(config.getQcStatusOutputPath(), geneToHaplotypeAnalysis);
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to create all output files", e);
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Peach");

        PeachConfig.addOptions(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PeachApplication peachApplication = new PeachApplication(configBuilder);
        peachApplication.run();
    }
}