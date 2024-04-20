package com.hartwig.hmftools.common.region;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class SpecificRegions
{
    public final List<ChrBaseRegion> Regions;
    public final List<String> Chromosomes;

    // config loading to filter on specific regions
    public static final String SPECIFIC_REGIONS = "specific_regions";
    public static final String SPECIFIC_REGIONS_DESC =
            "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd";

    public static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    public static final String SPECIFIC_CHROMOSOMES_DESC = "Restrict to chromosome(s) separated by ';'";

    private static final String SUB_ITEM_DELIM = ":";
    private static final String POS_ITEM_DELIM = "-";

    private static final Logger LOGGER = LogManager.getLogger(SpecificRegions.class);

    public SpecificRegions()
    {
        Chromosomes = Lists.newArrayList();
        Regions = Lists.newArrayList();
    }

    public boolean hasFilters() { return !Regions.isEmpty() || !Chromosomes.isEmpty(); }

    public void addRegion(final ChrBaseRegion newRegion)
    {
        if(!Chromosomes.contains(newRegion.Chromosome))
            Chromosomes.add(newRegion.Chromosome);

        for(ChrBaseRegion existingRegion : Regions)
        {
            if(existingRegion.overlaps(newRegion))
            {
                existingRegion.setStart(min(existingRegion.start(), newRegion.start()));
                existingRegion.setEnd(max(existingRegion.end(), newRegion.end()));
                return;
            }
        }

        Regions.add(newRegion);
    }

    public boolean includeChromosome(final String chromosome) { return Chromosomes.isEmpty() || Chromosomes.contains(chromosome); }
    public boolean excludeChromosome(final String chromosome) { return !includeChromosome(chromosome); }

    public boolean includeRegion(final ChrBaseRegion region)
    {
        return Regions.isEmpty() || Regions.stream().anyMatch(x -> x.overlaps(region));
    }

    public boolean includePosition(final String chromosome, final int position)
    {
        return Regions.isEmpty() || Regions.stream().anyMatch(x -> x.containsPosition(chromosome, position));
    }

    public boolean excludePosition(final String chromosome, final int position) { return !includePosition(chromosome, position); }

    public boolean includeRegion(final int posStart, final int posEnd)
    {
        // assumes chromosomes match
        return Regions.isEmpty() || Regions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), posStart, posEnd));
    }

    public boolean excludeRegion(final int posStart, final int posEnd) { return !includeRegion(posStart, posEnd); }

    public void log()
    {
        if(!Regions.isEmpty())
        {
            for(ChrBaseRegion region : Regions)
            {
                LOGGER.info("filtering for specific region: {}", region);
            }
        }
        else if(!Chromosomes.isEmpty())
        {
            LOGGER.info("filtering for specific chromosomes: {}", Chromosomes.stream().collect(Collectors.joining(",")));
        }
    }

    public static void addSpecificChromosomesRegionsConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SPECIFIC_CHROMOSOMES, SPECIFIC_CHROMOSOMES_DESC);
        configBuilder.addConfigItem(SPECIFIC_REGIONS, SPECIFIC_REGIONS_DESC);
    }

    public static SpecificRegions from(final ConfigBuilder configBuilder)
    {
        SpecificRegions specificRegions = new SpecificRegions();

        String regionsStr = configBuilder.getValue(SPECIFIC_REGIONS);
        String chromosomesStr = configBuilder.getValue(SPECIFIC_CHROMOSOMES);

        try
        {
            loadSpecificChromsomesOrRegions(regionsStr, chromosomesStr, specificRegions.Chromosomes, specificRegions.Regions);
            specificRegions.log();
        }
        catch(Exception e)
        {
            LOGGER.error("invalid specific config: regions({}) chromosomes({})", regionsStr, chromosomesStr);
            return null;
        }

        return specificRegions;
    }

    public static void loadSpecificChromsomesOrRegions(
            final ConfigBuilder configBuilder, final List<String> chromosomes, final List<ChrBaseRegion> regions) throws ParseException
    {
        loadSpecificChromsomesOrRegions(
                configBuilder.getValue(SPECIFIC_REGIONS), configBuilder.getValue(SPECIFIC_CHROMOSOMES), chromosomes, regions);
    }

    public static void loadSpecificChromsomesOrRegions(
            @Nullable final String specificRegionsStr, @Nullable final String specificChromosomesStr,
            final List<String> chromosomes, final List<ChrBaseRegion> regions) throws ParseException
    {
        if(specificRegionsStr != null && !specificRegionsStr.isEmpty())
        {
            loadSpecificRegions(specificRegionsStr, regions);

            Collections.sort(regions);

            for(ChrBaseRegion region : regions)
            {
                if(!chromosomes.contains(region.Chromosome))
                    chromosomes.add(region.Chromosome);
            }
        }
        else if(specificChromosomesStr != null && !specificChromosomesStr.isEmpty())
        {
            loadChromosomes(specificChromosomesStr, chromosomes);
        }
    }

    public static void loadChromosomes(final String chromosomesStr, final List<String> chromosomes)
    {
        if(chromosomesStr != null && !chromosomesStr.isEmpty())
        {
            Arrays.stream(chromosomesStr.split(ITEM_DELIM)).forEach(x -> chromosomes.add(x));
        }
    }

    public static List<String> loadSpecificChromsomes(final ConfigBuilder configBuilder)
    {
        return loadSpecificChromsomes(configBuilder.getValue(SPECIFIC_CHROMOSOMES));
    }

    public static List<String> loadSpecificChromsomes(final String specificChromosomesStr)
    {
        if(specificChromosomesStr == null || specificChromosomesStr.isEmpty())
            return Collections.emptyList();

        return Arrays.stream(specificChromosomesStr.split(ITEM_DELIM)).collect(Collectors.toList());
    }

    public static List<ChrBaseRegion> loadSpecificRegions(final ConfigBuilder configBuilder) throws ParseException
    {
        return loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS));
    }

    public static List<ChrBaseRegion> loadSpecificRegions(final String specificRegionsStr) throws ParseException
    {
        List<ChrBaseRegion> regions = Lists.newArrayList();
        loadSpecificRegions(specificRegionsStr, regions);
        return regions;
    }

    public static void loadSpecificRegions(final String specificRegionsStr, final List<ChrBaseRegion> regions) throws ParseException
    {
        if(specificRegionsStr == null || specificRegionsStr.isEmpty())
            return;

        final List<String> regionStrs = Arrays.stream(specificRegionsStr.split(ITEM_DELIM, -1)).collect(Collectors.toList());
        for(String regionStr : regionStrs)
        {
            ChrBaseRegion region = null;

            try
            {
                region = parseStandardFormat(regionStr);
            }
            catch(Exception e)
            {
                throw new ParseException(String.format("invalid specific region: %s", regionStr));
            }

            if(region == null || !region.isValid())
                throw new ParseException(String.format("invalid specific region: %s", region));

            regions.add(region);
        }
    }

    public static ChrBaseRegion parseStandardFormat(final String regionStr)
    {
        // expected format: chr:posStart-posEnd, separated by ';' or chr:posStart:posEnd
        final String[] items = regionStr.split(SUB_ITEM_DELIM);

        if(!(items.length == 3 || (items.length == 2 && items[1].contains(POS_ITEM_DELIM))))
        {
            return null;
        }

        String chromosome = items[0].intern();
        int posStart;
        int posEnd;

        if(items.length == 3)
        {
            posStart = Integer.parseInt(items[1]);
            posEnd = Integer.parseInt(items[2]);
        }
        else
        {
            String[] positions = items[1].split(POS_ITEM_DELIM, 2);
            posStart = Integer.parseInt(positions[0]);
            posEnd = Integer.parseInt(positions[1]);
        }

        return new ChrBaseRegion(chromosome, posStart, posEnd);
    }
}
