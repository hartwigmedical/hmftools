package com.hartwig.hmftools.cobalt.e2e;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.Nullable;

public abstract class SectionalData<T extends ChromosomeSection>
{
    private final SortedMap<HumanChromosome, List<T>> sectionsByChromosome = new TreeMap<>();

    @Nullable
    public abstract String header();

    public void addSection(T section)
    {
        List<T> sectionsForChromosome = sectionsByChromosome.computeIfAbsent(section.chromosome(), k -> new ArrayList<>());
        sectionsForChromosome.add(section);
    }

    public void write(File destination) throws IOException
    {
        List<String> lines = new ArrayList<>();
        String header = header();
        if (header != null)
        {
            lines.add(header);
        }
        sectionsByChromosome.forEach((chromosome, sections) -> sections.forEach(section -> lines.addAll(section.lines())));
        Files.write(destination.toPath(), lines);
    }
}

interface ChromosomeSection
{
    HumanChromosome chromosome();

    int start();

    int stop();

    default boolean is1Based()
    {
        return false;
    }

    String line(int position);

    default List<String> lines()
    {
        int interval = stop() - start();
        Preconditions.checkArgument(interval > 0);
        Preconditions.checkArgument(interval % 1000 == 0);

        int steps = (interval / 1000) + 1;
        List<java.lang.String> lines = new ArrayList<>();
        int position = start();
        for(int i = 0; i < steps; i++)
        {
            int bumpedPosition = is1Based() ? position + 1 : position;
            lines.add(line(bumpedPosition));
            position += 1000;
        }
        return lines;
    }
}
