package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;
import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

public class PeachApplication
{
    static final String CHAIN_FILE_DELIMITER = " ";
    static final String BED_FILE_DELIMITER = "\t";
    static final String TSV_DELIMITER = "\t";
    static final String HAPLOTYPE_EVENT_DELIMITER = ",";

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
        HaplotypePanel haplotypePanel = loadHaplotypePanel(config.haplotypesTsv);
        String callInputVcf;
        if (config.doLiftOver)
        {
            callInputVcf = getExtendedFileName(config.vcfFile, "liftover", ".vcf");
            String rejectVcf = getExtendedFileName(config.vcfFile, "reject", ".vcf");
            doLiftover(callInputVcf, rejectVcf);

            PCH_LOGGER.info("read bed of important regions");
            Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions = loadBedFile(config.liftOverBed);

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

            //TODO: validation of haplotype list making sense (try to do this in place where it fails early, but don't have to redo it during actual haplotype calling.)
        } else {
            callInputVcf = config.vcfFile;
        }
        Map<VariantHaplotypeEvent, Integer> variantEventToCount = loadRelevantVariantHaplotypeEvents(
                callInputVcf, haplotypePanel.getRelevantVariantPositions()
        );

        PCH_LOGGER.info("events found: {}", variantEventToCount.toString());

        Map<HaplotypeEvent, Set<String>> eventToRelevantGenes = variantEventToCount.keySet().stream()
                .collect(Collectors.toMap(e -> e, haplotypePanel::getRelevantGenes));

        PCH_LOGGER.info("genes per event: {}", eventToRelevantGenes.toString());

        PCH_LOGGER.info("finished running PEACH");
    }

    private Map<VariantHaplotypeEvent, Integer> loadRelevantVariantHaplotypeEvents(String vcf, Map<Chromosome, Set<Integer>> relevantVariantPositions)
    {
        Map<VariantHaplotypeEvent, Integer> eventToCount = new HashMap<>();
        try(
                AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                        vcf, new VCFCodec(), false)
        )
        {
            for(VariantContext variantContext : reader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;
                Chromosome chromosome = HumanChromosome.fromString(variantContext.getContig());
                if(!relevantVariantPositions.containsKey(chromosome))
                    continue;
                VariantHaplotypeEvent event = VariantHaplotypeEvent.fromVariantContext(variantContext);
                Set<Integer> relevantPositionsInChromosome = relevantVariantPositions.get(chromosome);
                if(event.getCoveredPositions().stream().noneMatch(relevantPositionsInChromosome::contains))
                    continue;
                Integer count = getEventCount(variantContext, event.id());
                if(count == 0)
                    continue;
                if(eventToCount.containsKey(event))
                {
                    PCH_LOGGER.error("encountered event with ID '{}' more than once in VCF '{}'", event.id(), vcf);
                    System.exit(1);
                }

                eventToCount.put(event, count);
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("failed to read VCF file({}): {}", vcf, e.toString());
            System.exit(1);
        }
        return eventToCount;
    }

    private Integer getEventCount(VariantContext variantContext, String eventId)
    {
        Integer count = null;
        GenotypeType genotype = variantContext.getGenotype(config.sampleName).getType();
        switch(genotype)
        {
            case NO_CALL:
            case HOM_REF:
                count = 0;
                break;
            case HET:
                count = 1;
                break;
            case HOM_VAR:
                count = 2;
                break;
            default:
                PCH_LOGGER.error("cannot get occurrence count for event with ID '{}' with genotype '{}'", eventId, genotype.toString());
                System.exit(1);
        }
        return count;
    }

    private void doLiftover(String liftOverVcf, String rejectVcf)
    {
        PCH_LOGGER.info("create adjusted chain file");
        String adjustedChainFile = getAdjustedChainFile(config.chainFile);

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

    private String getAdjustedChainFile(String chainFile)
    {
        String adjustedChainFile = getExtendedFileName(chainFile, "adjusted", ".over");
        try (
            Stream<String> lines = Files.lines(Paths.get(chainFile));
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

        return adjustedChainFile;
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

    private String getExtendedFileName(String originalFileName, String addition, String addBefore)
    {
        String[] fileItems = originalFileName.split("/");
        String filename = fileItems[fileItems.length - 1];
        int extensionIndex = filename.indexOf(addBefore);
        return config.outputDir + filename.substring(0, extensionIndex) + "." + addition + filename.substring(extensionIndex);
    }

    private HaplotypePanel loadHaplotypePanel(String filename)
    {
        List<Haplotype> haplotypes = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIMITER);

            int geneIndex = fieldsIndexMap.get("Gene");
            int haplotypeIndex = fieldsIndexMap.get("Haplotype");
            int wildTypeIndex = fieldsIndexMap.get("WildType");
            int eventsIndex = fieldsIndexMap.get("Events");

            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIMITER, -1);

                String haplotypeEventsString = values[eventsIndex];
                ImmutableSet<HaplotypeEvent> haplotypeEvents;
                if (haplotypeEventsString.isEmpty())
                    haplotypeEvents = ImmutableSet.of();
                else
                    haplotypeEvents = getHaplotypeEvents(haplotypeEventsString);

                Haplotype haplotype = new Haplotype(
                        values[geneIndex],
                        values[haplotypeIndex],
                        Boolean.parseBoolean(values[wildTypeIndex]),
                        haplotypeEvents
                );
                haplotypes.add(haplotype);
            }

            PCH_LOGGER.info("loaded {} haplotypes from file ({})", haplotypes.size(), filename);
        }
        catch(Exception e)
        {
            PCH_LOGGER.error("failed to load haplotypes TSV({}): {}", filename, e.toString());
            System.exit(1);
        }
        return new HaplotypePanel(haplotypes);
    }

    private static ImmutableSet<HaplotypeEvent> getHaplotypeEvents(String haplotypeEventsString)
    {
        return Arrays.stream(haplotypeEventsString.split(HAPLOTYPE_EVENT_DELIMITER))
                .map(HaplotypeEventFactory::fromId).collect(ImmutableSet.toImmutableSet());
    }

    private Map<Chromosome, List<BaseRegion>> loadBedFile(final String filename)
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