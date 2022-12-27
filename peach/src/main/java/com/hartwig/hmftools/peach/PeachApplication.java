package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.peach.data_loader.HaplotypeEventLoader;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.output.AllHaplotypeCombinationsFile;
import com.hartwig.hmftools.peach.output.EventsFile;
import com.hartwig.hmftools.peach.output.EventsPerGeneFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Stream;

import static com.hartwig.hmftools.peach.PeachUtils.BED_FILE_DELIMITER;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;
import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

public class PeachApplication
{
    public static final String CHAIN_FILE_DELIMITER = " ";

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
            callInputVcf = config.getLiftoverOutputVcfPath();
            String rejectVcf = config.getLiftoverRejectVcfPath();
            doLiftover(callInputVcf, rejectVcf);

            PCH_LOGGER.info("read bed of important regions");
            Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions = loadRegionsToLiftover(config.liftOverBed);

            PCH_LOGGER.info("check rejected liftover variants for relevance");
            int potentiallyMissedCount = countPotentiallyRelevantVariantsMissed(rejectVcf, chromosomeToRelevantRegions);

            if (potentiallyMissedCount == 0)
            {
                PCH_LOGGER.info("all potentially relevant variants have been lifted over");
            }
            else
            {
                PCH_LOGGER.warn("some potentially relevant variants have not been lifted over: {}", potentiallyMissedCount);
            }
            //TODO: handle unlifted relevant variants properly
            //TODO: handle reference sequence differences V37 vs V38 properly
            //TODO: maybe force at most two haplotype calls per gene in output. Maybe do this optionally. Maybe by default
        } else {
            callInputVcf = config.vcfFile;
        }
        Map<String, Integer> eventIdToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                callInputVcf, config.sampleName, haplotypePanel.getRelevantVariantPositions()
        );

        PCH_LOGGER.info("events found: {}", eventIdToCount.toString());

        HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = caller.getGeneToHaplotypeAnalysis(eventIdToCount);

        PCH_LOGGER.info("haplotypes called: {}", geneToHaplotypeAnalysis.toString());

        PCH_LOGGER.info("Write events output file");
        try
        {
            EventsFile.write(config.getEventsOutputPath(), eventIdToCount);
        }
        catch (IOException e)
        {
            PCH_LOGGER.error("failed to create events output file({}): {}", config.getEventsOutputPath(), e.toString());
            System.exit(1);
        }

        PCH_LOGGER.info("Write events per gene output file");
        try
        {
            EventsPerGeneFile.write(config.getEventsPerGeneOutputPath(), geneToHaplotypeAnalysis);
        }
        catch (IOException e)
        {
            PCH_LOGGER.error("failed to create events per gene output file({}): {}", config.getEventsPerGeneOutputPath(), e.toString());
            System.exit(1);
        }
        PCH_LOGGER.info("Write all haplotype combinations output file");
        try
        {
            AllHaplotypeCombinationsFile.write(config.getAllHaplotypeCombinationsOutputPath(), geneToHaplotypeAnalysis);
        }
        catch (IOException e)
        {
            PCH_LOGGER.error("failed to create all haplotype combinations output file({}): {}", config.getAllHaplotypeCombinationsOutputPath(), e.toString());
            System.exit(1);
        }

        PCH_LOGGER.info("finished running PEACH");
    }

    private void doLiftover(String liftOverVcf, String rejectVcf)
    {
        PCH_LOGGER.info("create adjusted chain file");
        String adjustedChainFile = config.getAdjustedChainFilePath();
        adjustChainFile(config.chainFile, adjustedChainFile);

        PCH_LOGGER.info("do lift over");
        ProcessBuilder pb = new ProcessBuilder(
                "java",
                "-jar",
                config.picardJar,
                "LiftoverVcf",
                "CHAIN=" + adjustedChainFile,
                "INPUT=" + config.vcfFile,
                "OUTPUT=" + liftOverVcf,
                "REFERENCE_SEQUENCE=" + config.targetRefGenome,
                "REJECT=" + rejectVcf,
                "RECOVER_SWAPPED_REF_ALT=true",
                "WRITE_ORIGINAL_POSITION=true",
                "WRITE_ORIGINAL_ALLELES=true"
        );
        try
        {
            pb.inheritIO();
            Process process = pb.start();
            int exitCode = process.waitFor();
            if (exitCode != 0)
            {
                PCH_LOGGER.error("Picard had a non-zero exit code: {}", exitCode);
                System.exit(1);
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("Picard LiftoverVcf failed: ");
            e.printStackTrace();
            System.exit(1);
        }
        catch(InterruptedException e)
        {
            PCH_LOGGER.error("Picard LiftoverVcf was interrupted");
            e.printStackTrace();
            System.exit(1);
        }
    }

    private int countPotentiallyRelevantVariantsMissed(String rejectVcf, Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions)
    {
        int potentiallyMissedVariantCount = 0;
        try(
                AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                        rejectVcf, new VCFCodec(), false)
        )
        {
            for(VariantContext variantContext : reader.iterator())
            {
                if (isPotentiallyRelevant(variantContext, chromosomeToRelevantRegions))
                {
                    potentiallyMissedVariantCount += 1;
                }
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("failed to read rejected liftover VCF file({}): {}", rejectVcf, e.toString());
            System.exit(1);
        }
        return potentiallyMissedVariantCount;
    }

    private boolean isPotentiallyRelevant(VariantContext variantContext, Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions)
    {
        Chromosome variantChromosome = HumanChromosome.fromString(variantContext.getContig());
        if (!chromosomeToRelevantRegions.containsKey(variantChromosome))
            return false;

        BaseRegion variantRegion = new BaseRegion(variantContext.getStart(), variantContext.getEnd());

        return chromosomeToRelevantRegions.get(variantChromosome).stream().anyMatch(r -> r.overlaps(variantRegion));
    }

    private void adjustChainFile(String sourceChainFile, String adjustedChainFile)
    {
        try (
            Stream<String> lines = Files.lines(Paths.get(sourceChainFile));
            PrintWriter pw = new PrintWriter(adjustedChainFile, StandardCharsets.UTF_8)
        )
        {
            lines.forEachOrdered(line-> pw.println(getAdjustedChainFileLine(line)));
        }
        catch (IOException e)
        {
            PCH_LOGGER.error("could not create adjusted chain file: ");
            e.printStackTrace();
            System.exit(1);
        }
    }

    private Map<Chromosome, List<BaseRegion>> loadRegionsToLiftover(final String filename)
    {
        final Map<Chromosome,List<BaseRegion>> chromosomeToRegions = Maps.newHashMap();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            for(String line : lines)
            {
                final String[] values = line.split(BED_FILE_DELIMITER, -1);

                Chromosome chromosome = HumanChromosome.fromString(values[0]);
                int posStart = Integer.parseInt(values[1]) + 1; // as per convention
                int posEnd = Integer.parseInt(values[2]);

                if(!chromosomeToRegions.containsKey(chromosome))
                {
                    chromosomeToRegions.put(chromosome, Lists.newArrayList());
                }

                chromosomeToRegions.get(chromosome).add(new BaseRegion(posStart, posEnd));
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("failed to load BED file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return chromosomeToRegions;
    }

    private String getAdjustedChainFileLine(String line)
    {
        if (line.startsWith("chain"))
        {
            String[] items = line.split(CHAIN_FILE_DELIMITER);
            StringJoiner newLineJoiner = new StringJoiner(CHAIN_FILE_DELIMITER);
            for (int i = 0; i < items.length; i++)
            {
                if (i == 2)
                    newLineJoiner.add(RefGenomeFunctions.stripChrPrefix(items[i]));
                else
                    newLineJoiner.add(items[i]);
            }
            return newLineJoiner.toString();
        }
        else
            return line;
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