package com.hartwig.hmftools.peach;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Streams;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.NonWildTypeHaplotype;
import com.hartwig.hmftools.peach.haplotype.WildTypeHaplotype;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;
import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

public class DataLoader
{
    public static final String TSV_DELIMITER = "\t";
    public static final String HAPLOTYPE_EVENT_DELIMITER = ",";
    public static final String BED_FILE_DELIMITER = "\t";

    public static Map<String, Integer> loadRelevantVariantHaplotypeEvents(
            String vcf, String sampleName, Map<Chromosome, Set<Integer>> relevantVariantPositions
    )
    {
        Map<String, Integer> eventIdToCount = new HashMap<>();
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
                Integer count = getEventCount(variantContext.getGenotype(sampleName).getType(), event.id());
                if(count == 0)
                    continue;

                if(eventIdToCount.containsKey(event.id()))
                {
                    PCH_LOGGER.error("encountered event with ID '{}' more than once in VCF '{}'", event.id(), vcf);
                    System.exit(1);
                }

                eventIdToCount.put(event.id(), count);
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("failed to read VCF file({}): {}", vcf, e.toString());
            System.exit(1);
        }
        return eventIdToCount;
    }

    public static HaplotypePanel loadHaplotypePanel(String filename)
    {
        Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel = loadGeneToGeneHaplotypePanel(filename);
        return new HaplotypePanel(geneToGeneHaplotypePanel);
    }

    public static Map<Chromosome, List<BaseRegion>> loadBedFile(final String filename)
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

    private static Integer getEventCount(GenotypeType genotypeType, String eventId)
    {
        Integer count = null;
        switch(genotypeType)
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
                PCH_LOGGER.error("cannot get occurrence count for event with ID '{}' with genotypeType '{}'", eventId, genotypeType.toString());
                System.exit(1);
        }
        return count;
    }

    @NotNull
    private static Map<String, GeneHaplotypePanel> loadGeneToGeneHaplotypePanel(String filename)
    {
        Map<String, List<WildTypeHaplotype>> geneToWildTypeHaplotypes = new HashMap<>();
        Map<String, List<NonWildTypeHaplotype>> geneToNonWildTypeHaplotypes = new HashMap<>();
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
                String gene = values[geneIndex];
                boolean isWildType = Boolean.parseBoolean(values[wildTypeIndex]);

                ImmutableList<HaplotypeEvent> haplotypeEvents = getHaplotypeEvents(haplotypeEventsString);
                if (isWildType)
                {
                    WildTypeHaplotype haplotype = new WildTypeHaplotype(values[haplotypeIndex], haplotypeEvents);
                    if (!geneToWildTypeHaplotypes.containsKey(gene))
                        geneToWildTypeHaplotypes.put(gene, Lists.newArrayList());
                    geneToWildTypeHaplotypes.get(gene).add(haplotype);
                }
                else
                {
                    NonWildTypeHaplotype haplotype = new NonWildTypeHaplotype(values[haplotypeIndex], haplotypeEvents);
                    if (!geneToNonWildTypeHaplotypes.containsKey(gene))
                        geneToNonWildTypeHaplotypes.put(gene, Lists.newArrayList());
                    geneToNonWildTypeHaplotypes.get(gene).add(haplotype);
                }
            }
        }
        catch(Exception e)
        {
            PCH_LOGGER.error("failed to load haplotypes TSV({}): {}", filename, e.toString());
            System.exit(1);
        }

        int haplotypeCount = Streams.concat(
                geneToWildTypeHaplotypes.values().stream(),
                geneToNonWildTypeHaplotypes.values().stream()
        ).mapToInt(List::size).sum();
        PCH_LOGGER.info("loaded {} haplotypes from file ({})", haplotypeCount, filename);

        return createGeneToGeneHaplotypePanel(geneToWildTypeHaplotypes, geneToNonWildTypeHaplotypes);
    }

    private static Map<String, GeneHaplotypePanel> createGeneToGeneHaplotypePanel(
            Map<String, List<WildTypeHaplotype>> geneToWildTypeHaplotypes,
            Map<String, List<NonWildTypeHaplotype>> geneToNonWildTypeHaplotypes
    )
    {
        Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel = new HashMap<>();
        Set<String> genes = Streams.concat(
                geneToWildTypeHaplotypes.keySet().stream(),
                geneToNonWildTypeHaplotypes.keySet().stream()
        ).collect(Collectors.toSet());
        for (String gene : genes)
        {
            List<WildTypeHaplotype> wildTypeHaplotypes = geneToWildTypeHaplotypes.getOrDefault(gene, new ArrayList<>());
            if (wildTypeHaplotypes.size() != 1)
            {
                throw new RuntimeException(String.format("Cannot have more than 1 wild type haplotype for gene: %s", gene));
            }
            WildTypeHaplotype wildTypeHaplotype = wildTypeHaplotypes.get(0);

            GeneHaplotypePanel geneHaplotypePanel = new GeneHaplotypePanel(
                    wildTypeHaplotype,
                    ImmutableList.copyOf(geneToNonWildTypeHaplotypes.getOrDefault(gene, new ArrayList<>()))
            );
            geneToGeneHaplotypePanel.put(gene, geneHaplotypePanel);
        }

        return geneToGeneHaplotypePanel;
    }

    private static ImmutableList<HaplotypeEvent> getHaplotypeEvents(String haplotypeEventsString)
    {
        if (haplotypeEventsString.isEmpty())
            return ImmutableList.of();
        else
            return Arrays.stream(haplotypeEventsString.split(HAPLOTYPE_EVENT_DELIMITER))
                    .map(HaplotypeEventFactory::fromId).collect(ImmutableList.toImmutableList());
    }
}
