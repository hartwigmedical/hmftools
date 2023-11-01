package com.hartwig.hmftools.peach;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.region.BaseRegion;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

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

public class LiftoverService
{
    public static final String CHAIN_FILE_DELIMITER = " ";

    private final PeachConfig config;

    public LiftoverService(PeachConfig config)
    {
        this.config = config;
    }

    public void doLiftover()
    {
        PCH_LOGGER.info("create adjusted chain file");
        adjustChainFile();

        PCH_LOGGER.info("run Picard LiftoverVcf");
        runPicardLiftover();

        PCH_LOGGER.info("check rejected liftover variants for relevance");
        int potentiallyMissedCount = countPotentiallyRelevantVariantsMissed();

        if (potentiallyMissedCount == 0)
        {
            PCH_LOGGER.info("all potentially relevant variants have been lifted over");
        }
        else
        {
            PCH_LOGGER.error("some potentially relevant variants have not been lifted over: {}", potentiallyMissedCount);
            System.exit(1);
        }
    }

    private void adjustChainFile()
    {
        try (
                Stream<String> lines = Files.lines(Paths.get(config.chainFile));
                PrintWriter pw = new PrintWriter(config.getAdjustedChainFilePath(), StandardCharsets.UTF_8)
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

    private void runPicardLiftover()
    {
        ProcessBuilder pb = new ProcessBuilder(
                "java",
                "-jar",
                config.picardJar,
                "LiftoverVcf",
                "CHAIN=" + config.getAdjustedChainFilePath(),
                "INPUT=" + config.vcfFile,
                "OUTPUT=" + config.getLiftoverOutputVcfPath(),
                "REFERENCE_SEQUENCE=" + config.targetRefGenome,
                "REJECT=" + config.getLiftoverRejectVcfPath(),
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

    private int countPotentiallyRelevantVariantsMissed()
    {
        Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions = loadRegionsToLiftover();
        int potentiallyMissedVariantCount = 0;
        try(
                AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                        config.getLiftoverRejectVcfPath(), new VCFCodec(), false)
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
            PCH_LOGGER.error("failed to read rejected liftover VCF file({}): {}", config.getLiftoverRejectVcfPath(), e.toString());
            System.exit(1);
        }
        return potentiallyMissedVariantCount;
    }

    private Map<Chromosome, List<BaseRegion>> loadRegionsToLiftover()
    {
        final Map<Chromosome,List<BaseRegion>> chromosomeToRegions = Maps.newHashMap();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(config.liftOverBed));
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
            PCH_LOGGER.error("failed to load BED file({}): {}", config.liftOverBed, e.toString());
            System.exit(1);
        }

        return chromosomeToRegions;
    }

    private boolean isPotentiallyRelevant(VariantContext variantContext, Map<Chromosome, List<BaseRegion>> chromosomeToRelevantRegions)
    {
        Chromosome variantChromosome = HumanChromosome.fromString(variantContext.getContig());
        if (!chromosomeToRelevantRegions.containsKey(variantChromosome))
            return false;

        BaseRegion variantRegion = new BaseRegion(variantContext.getStart(), variantContext.getEnd());

        return chromosomeToRelevantRegions.get(variantChromosome).stream().anyMatch(r -> r.overlaps(variantRegion));
    }
}
