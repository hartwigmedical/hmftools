package com.hartwig.hmftools.peach.data_loader;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.hartwig.hmftools.peach.panel.GeneHaplotypePanel;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;

import org.jetbrains.annotations.NotNull;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PanelLoader
{
    public static final String HAPLOTYPE_EVENT_DELIMITER = ";";

    @NotNull
    public static HaplotypePanel loadHaplotypePanel(@NotNull String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            HaplotypePanel panel = loadHaplotypePanel(lines);

            PCH_LOGGER.info("loaded {} haplotypes from file ({})", panel.getHaplotypeCount(), filename);
            return panel;
        }
        catch(Exception e)
        {
            throw new RuntimeException(String.format("failed to load haplotypes TSV: %s", filename), e);
        }
    }

    @NotNull
    private static HaplotypePanel loadHaplotypePanel(@NotNull List<String> lines)
    {
        Map<String, List<DefaultHaplotype>> geneToDefaultHaplotypes = new HashMap<>();
        Map<String, List<NonDefaultHaplotype>> geneToNonDefaultHaplotypes = new HashMap<>();

        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        int geneIndex = fieldsIndexMap.get("gene");
        int haplotypeIndex = fieldsIndexMap.get("haplotype");
        int defaultIndex = fieldsIndexMap.get("default");
        int wildTypeIndex = fieldsIndexMap.get("wildType");
        int eventsIndex = fieldsIndexMap.get("events");

        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(TSV_DELIM, -1);

            String haplotypeEventsString = values[eventsIndex];
            String gene = values[geneIndex];
            boolean isDefault = Boolean.parseBoolean(values[defaultIndex]);
            boolean isWildType = Boolean.parseBoolean(values[wildTypeIndex]);

            ImmutableList<HaplotypeEvent> haplotypeEvents = getHaplotypeEvents(haplotypeEventsString);
            if(isDefault)
            {
                DefaultHaplotype haplotype = new DefaultHaplotype(values[haplotypeIndex], isWildType, haplotypeEvents);
                if(!geneToDefaultHaplotypes.containsKey(gene))
                {
                    geneToDefaultHaplotypes.put(gene, Lists.newArrayList());
                }
                geneToDefaultHaplotypes.get(gene).add(haplotype);
            }
            else
            {
                NonDefaultHaplotype haplotype = new NonDefaultHaplotype(values[haplotypeIndex], isWildType, haplotypeEvents);
                if(!geneToNonDefaultHaplotypes.containsKey(gene))
                {
                    geneToNonDefaultHaplotypes.put(gene, Lists.newArrayList());
                }
                geneToNonDefaultHaplotypes.get(gene).add(haplotype);
            }
        }

        Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel =
                createGeneToGeneHaplotypePanel(geneToDefaultHaplotypes, geneToNonDefaultHaplotypes);
        return new HaplotypePanel(geneToGeneHaplotypePanel);
    }

    @NotNull
    private static Map<String, GeneHaplotypePanel> createGeneToGeneHaplotypePanel(
            @NotNull Map<String, List<DefaultHaplotype>> geneToDefaultHaplotypes,
            @NotNull Map<String, List<NonDefaultHaplotype>> geneToNonDefaultHaplotypes)
    {
        Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel = new HashMap<>();
        Set<String> genes = Streams.concat(geneToDefaultHaplotypes.keySet().stream(), geneToNonDefaultHaplotypes.keySet().stream())
                .collect(Collectors.toSet());
        for(String gene : genes)
        {
            List<DefaultHaplotype> defaultHaplotypes = geneToDefaultHaplotypes.getOrDefault(gene, new ArrayList<>());
            if(defaultHaplotypes.size() != 1)
            {
                throw new RuntimeException(String.format("Need to have exactly 1 default haplotype for gene: %s", gene));
            }
            DefaultHaplotype defaultHaplotype = defaultHaplotypes.get(0);
            List<NonDefaultHaplotype> nonDefaultHaplotypes = geneToNonDefaultHaplotypes.getOrDefault(gene, new ArrayList<>());
            if(nonDefaultHaplotypes.stream().anyMatch(h -> h.getName().equals(defaultHaplotype.getName())))
            {
                throw new RuntimeException(String.format("Cannot have non-default haplotype with same name as default haplotype for gene: %s", gene));
            }

            String wildTypeHaplotypeName = getWildTypeHaplotypeName(defaultHaplotype, nonDefaultHaplotypes, gene);

            GeneHaplotypePanel geneHaplotypePanel =
                    new GeneHaplotypePanel(defaultHaplotype, ImmutableList.copyOf(nonDefaultHaplotypes), wildTypeHaplotypeName);
            geneToGeneHaplotypePanel.put(gene, geneHaplotypePanel);
        }

        return geneToGeneHaplotypePanel;
    }

    @NotNull
    private static String getWildTypeHaplotypeName(@NotNull DefaultHaplotype defaultHaplotype,
            @NotNull List<NonDefaultHaplotype> nonDefaultHaplotypes, @NotNull String gene)
    {
        List<String> wildTypeHaplotypeNames = nonDefaultHaplotypes.stream()
                .filter(NonDefaultHaplotype::isWildType)
                .map(NonDefaultHaplotype::getName)
                .collect(Collectors.toList());
        if(defaultHaplotype.isWildType())
        {
            wildTypeHaplotypeNames.add(defaultHaplotype.getName());
        }

        if(wildTypeHaplotypeNames.size() != 1)
        {
            throw new RuntimeException(String.format("Need to have exactly 1 wild type haplotype for gene: %s", gene));
        }

        String wildTypeHaplotypeName = wildTypeHaplotypeNames.get(0);

        boolean haplotypeStatusesClash = !defaultHaplotype.isWildType() && nonDefaultHaplotypes.stream()
                .anyMatch(h -> h.getName().equals(wildTypeHaplotypeName) && !h.isWildType());
        if(haplotypeStatusesClash)
        {
            throw new RuntimeException(String.format("Wild type status inconsistent for haplotype %s for gene %s", wildTypeHaplotypeName, gene));
        }

        return wildTypeHaplotypeName;
    }

    @NotNull
    private static ImmutableList<HaplotypeEvent> getHaplotypeEvents(@NotNull String haplotypeEventsString)
    {
        if(haplotypeEventsString.isEmpty())
        {
            return ImmutableList.of();
        }
        else
        {
            return Arrays.stream(haplotypeEventsString.split(HAPLOTYPE_EVENT_DELIMITER))
                    .map(HaplotypeEventFactory::fromId)
                    .collect(ImmutableList.toImmutableList());
        }
    }
}
