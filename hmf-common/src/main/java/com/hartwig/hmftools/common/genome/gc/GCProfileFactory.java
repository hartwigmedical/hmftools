package com.hartwig.hmftools.common.genome.gc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.Options;

public final class GCProfileFactory
{
    private static final String RATIO_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int GC_CONTENT_COLUMN = 2;
    private static final int NON_N_PERCENTAGE_COLUMN = 3;
    private static final int MAPPABLE_PERCENTAGE_COLUMN = 4;

    // the file entries define the window size so a constant ought not be required
    public static final int WINDOW_SIZE = 1000;

    public static final String GC_PROFILE = "gc_profile";
    public static final String GC_PROFILE_DESC = "Path to GC profile";

    public static void addGcProfilePath(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPath(GC_PROFILE, required, GC_PROFILE_DESC);
    }

    public static Multimap<Chromosome, GCProfile> loadGCContent(final String fileName) throws IOException
    {
        return loadGCContent(WINDOW_SIZE, Files.readAllLines(new File(fileName).toPath()));
    }

    public static Map<String,List<GCProfile>> loadChrGcProfileMap(final String fileName) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(fileName).toPath());
        Map<String,List<GCProfile>> chrProfileMap = Maps.newHashMap();

        String currentChr = "";
        List<GCProfile> profiles = null;

        for(String line : lines)
        {
            final GCProfile gcProfile = fromLine(WINDOW_SIZE, line);
            if(!gcProfile.chromosome().equals(currentChr))
            {
                profiles = Lists.newArrayList();
                currentChr = gcProfile.chromosome();
                chrProfileMap.put(gcProfile.chromosome(), profiles);
            }

            chrProfileMap.put(gcProfile.chromosome(), profiles);
            profiles.add(gcProfile);
        }

        return chrProfileMap;
    }

    public static ListMultimap<Chromosome, GCProfile> loadGCContent(int windowSize, final String fileName) throws IOException
    {
        return loadGCContent(windowSize, Files.readAllLines(new File(fileName).toPath()));
    }

    private static ListMultimap<Chromosome, GCProfile> loadGCContent(int windowSize, final List<String> lines)
    {
        final ListMultimap<Chromosome, GCProfile> result = ArrayListMultimap.create();
        for(String line : lines)
        {
            final GCProfile gcProfile = fromLine(windowSize, line);
            if(HumanChromosome.contains(gcProfile.chromosome()))
            {
                result.put(HumanChromosome.fromString(gcProfile.chromosome()), gcProfile);
            }
        }

        return result;
    }

    private static GCProfile fromLine(int windowSize, final String ratioLine)
    {
        final String[] values = ratioLine.split(RATIO_COLUMN_SEPARATOR);

        String chromosome = values[CHROMOSOME_COLUMN].trim();
        int position = Integer.parseInt(values[START_FIELD_COLUMN].trim());
        double gcContent = Double.parseDouble(values[GC_CONTENT_COLUMN].trim());
        double nonNPercentage = Double.parseDouble(values[NON_N_PERCENTAGE_COLUMN].trim());
        double mappablePercentage = Double.parseDouble(values[MAPPABLE_PERCENTAGE_COLUMN].trim());

        return ImmutableGCProfile.builder()
                .chromosome(chromosome)
                .start(position + 1) // GCProfile is zero-indexed
                .end(position + windowSize)
                .gcContent(gcContent)
                .nonNPercentage(nonNPercentage)
                .mappablePercentage(mappablePercentage)
                .build();
    }
}
