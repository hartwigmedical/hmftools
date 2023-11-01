package com.hartwig.hmftools.peach.data_loader;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.hartwig.hmftools.peach.panel.GeneHaplotypePanel;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;
import com.hartwig.hmftools.peach.PeachUtils;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import com.hartwig.hmftools.peach.haplotype.NonWildTypeHaplotype;
import com.hartwig.hmftools.peach.haplotype.WildTypeHaplotype;
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

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class PanelLoader
{
    public static HaplotypePanel loadHaplotypePanel(String filename)
    {
        HaplotypePanel panel = null;
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            panel = loadHaplotypePanel(lines);

            PCH_LOGGER.info("loaded {} haplotypes from file ({})", panel.getHaplotypeCount(), filename);
        }
        catch(Exception e)
        {
            PCH_LOGGER.error("failed to load haplotypes TSV({}): {}", filename, e.toString());
            System.exit(1);
        }
        return panel;
    }

    @NotNull
    private static HaplotypePanel loadHaplotypePanel(List<String> lines)
    {
        Map<String, List<WildTypeHaplotype>> geneToWildTypeHaplotypes = new HashMap<>();
        Map<String, List<NonWildTypeHaplotype>> geneToNonWildTypeHaplotypes = new HashMap<>();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, PeachUtils.TSV_DELIMITER);
        int geneIndex = fieldsIndexMap.get("Gene");
        int haplotypeIndex = fieldsIndexMap.get("Haplotype");
        int wildTypeIndex = fieldsIndexMap.get("WildType");
        int eventsIndex = fieldsIndexMap.get("Events");

        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(PeachUtils.TSV_DELIMITER, -1);

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

        Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel = createGeneToGeneHaplotypePanel(
                geneToWildTypeHaplotypes, geneToNonWildTypeHaplotypes
        );
        return new HaplotypePanel(geneToGeneHaplotypePanel);
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
            return Arrays.stream(haplotypeEventsString.split(PeachUtils.HAPLOTYPE_EVENT_DELIMITER))
                    .map(HaplotypeEventFactory::fromId).collect(ImmutableList.toImmutableList());
    }
}
