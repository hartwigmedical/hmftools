package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import com.hartwig.hmftools.peach.data_loader.HaplotypeEventLoader;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.output.AllHaplotypeCombinationsFile;
import com.hartwig.hmftools.peach.output.BestHaplotypeCombinationsFile;
import com.hartwig.hmftools.peach.output.EventsFile;
import com.hartwig.hmftools.peach.output.EventsPerGeneFile;
import com.hartwig.hmftools.peach.output.QcStatusFile;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PeachApplication
{
    @NotNull
    private final PeachConfig config;

    public PeachApplication(@NotNull final PeachConfig config)
    {
        this.config = config;
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
        final Options options = PeachConfig.createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            PeachConfig config = new PeachConfig(cmd);
            PeachApplication peachApplication = new PeachApplication(config);
            peachApplication.run();
        }
        catch(ParseException e)
        {
            PCH_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PeachApplication", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}