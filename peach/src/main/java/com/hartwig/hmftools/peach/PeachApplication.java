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
        if(!config.isValid())
        {
            PCH_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }
        PCH_LOGGER.info("running PEACH");

        PCH_LOGGER.info("creating output directory");
        // Output directory cannot be null in valid config
        assert config.outputDir != null;
        File outputDirectory = new File(config.outputDir);
        if (!outputDirectory.exists() && !outputDirectory.mkdirs())
        {
            PCH_LOGGER.error("could not create output directory, exiting");
            System.exit(1);
        }

        PCH_LOGGER.info("read haplotypes TSV");
        HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(config.haplotypesTsv);
        String callInputVcf;
        if (config.doLiftOver)
        {

            LiftoverService liftoverService = new LiftoverService(config);
            liftoverService.doLiftover();
            callInputVcf = config.getLiftoverOutputVcfPath();
            //TODO: handle reference sequence differences V37 vs V38 properly
        } else {
            callInputVcf = config.vcfFile;
        }
        PCH_LOGGER.info("read events");
        Map<String, Integer> eventIdToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                callInputVcf, config.sampleName, haplotypePanel.getRelevantVariantPositions()
        );

        PCH_LOGGER.info("call haplotypes");
        HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = caller.getGeneToHaplotypeAnalysis(eventIdToCount);

        writeOutputFiles(eventIdToCount, geneToHaplotypeAnalysis);

        PCH_LOGGER.info("finished running PEACH");
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
        catch (IOException e)
        {
            PCH_LOGGER.error("failed to create all output files: {}", e.toString());
            System.exit(1);
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