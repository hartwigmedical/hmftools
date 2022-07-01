package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ChrBaseRegion implements Cloneable, Comparable<ChrBaseRegion>
{
    public final String Chromosome;
    public int[] Positions;

    public ChrBaseRegion(final String chromosome, final int[] positions)
    {
        Chromosome = chromosome;
        Positions = positions;
    }

    public ChrBaseRegion(final String chromosome, final int posStart, final int posEnd)
    {
        Chromosome = chromosome;
        Positions = new int[] { posStart, posEnd };
    }

    public static ChrBaseRegion from(final GenomeRegion region) { return new ChrBaseRegion(region.chromosome(), region.start(), region.end()); }

    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }
    public String chromosome() { return Chromosome; }

    public void setPosition(int position, int index) { Positions[index] = position; }
    public void setStart(int pos) { setPosition(pos, SE_START); }
    public void setEnd(int pos) { setPosition(pos, SE_END); }

    public int baseLength() { return length() + 1; }
    public int length() { return Positions[SE_END] - Positions[SE_START]; }

    public boolean isValid() { return HumanChromosome.contains(Chromosome) && hasValidPositions(); }
    public boolean hasValidPositions() { return Positions[SE_START] > 0 & Positions[SE_END] >= Positions[SE_START]; }

    public boolean overlaps(final ChrBaseRegion other)
    {
        if(!Chromosome.equals(other.Chromosome))
            return false;

        return positionsOverlap(Positions[SE_START], Positions[SE_END], other.Positions[SE_START], other.Positions[SE_END]);
    }

    public boolean containsPosition(int position) { return positionWithin(position, start(), end()); }

    public boolean containsPosition(final String chromosome, int position)
    {
        return Chromosome.equals(chromosome) && positionWithin(position, start(), end());
    }

    public boolean matches(final ChrBaseRegion other)
    {
        return Chromosome.equals(other.Chromosome) && start() == other.start() && end() == other.end();
    }

    public String toString() { return String.format("%s:%d-%d", Chromosome, Positions[SE_START], Positions[SE_END]); }

    @Override
    public Object clone()
    {
        try
        {
            ChrBaseRegion br = (ChrBaseRegion) super.clone();
            br.Positions = Positions.clone();
            return br;
        }
        catch (CloneNotSupportedException e)
        {
            // Will not happen in this case
            return null;
        }
    }

    @Override
    public boolean equals(Object obj)
    {
        // same instance
        if (obj == this) { return true; }
        // null
        if (obj == null) { return false; }
        // type
        if (!getClass().equals(obj.getClass())) { return false; }
        // cast and compare state
        ChrBaseRegion other = (ChrBaseRegion) obj;
        return matches(other);
    }

    @Override
    public int compareTo(@NotNull final ChrBaseRegion other)
    {
        if(Chromosome.equals(other.Chromosome))
        {
            if (start() < other.start())
            {
                return -1;
            }
            else if (start() == other.start())
            {
                return 0;
            }
            return 1;
        }

        return ContigComparator.INSTANCE.compare(Chromosome, other.Chromosome);
    }

    // config loading to filter on specific regions
    public static final String SPECIFIC_REGIONS = "specific_regions";
    public static final String SPECIFIC_REGIONS_DESC =
            "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd";

    public static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    public static final String SPECIFIC_CHROMOSOMES_DESC = "Restrict to chromosome(s) separated by ';'";

    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = ":";
    public static final String POS_ITEM_DELIM = "-";

    public static void addSpecificChromosomesRegionsConfig(final Options options)
    {
        options.addOption(SPECIFIC_CHROMOSOMES, true, SPECIFIC_CHROMOSOMES_DESC);
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
    }

    public static List<ChrBaseRegion> loadSpecificRegions(final CommandLine cmd) throws ParseException
    {
        List<ChrBaseRegion> regions = Lists.newArrayList();

        if(!cmd.hasOption(SPECIFIC_REGIONS))
            return regions;

        // expected format: chromosome:positionStart:positionEnd, separated by ';'
        final List<String> regionStrs = Arrays.stream(cmd.getOptionValue(SPECIFIC_REGIONS).split(ITEM_DELIM, -1)).collect(Collectors.toList());
        for(String regionStr : regionStrs)
        {
            final String[] items = regionStr.split(SUB_ITEM_DELIM);

            if(!(items.length == 3 || (items.length == 2 && items[1].contains(POS_ITEM_DELIM))))
            {
                throw new ParseException(String.format("invalid specific region: %s", regionStr));
            }

            String chromosome = items[0];
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

            ChrBaseRegion region = new ChrBaseRegion(chromosome, posStart, posEnd);

            if(!region.isValid())
                throw new ParseException(String.format("invalid specific region: %s", region));

            regions.add(region);
        }

        return regions;
    }

    public static List<String> loadSpecificChromsomes(final CommandLine cmd)
    {
        return cmd.hasOption(SPECIFIC_CHROMOSOMES) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_CHROMOSOMES).split(ITEM_DELIM)).collect(Collectors.toList()) : Lists.newArrayList();
    }

    public static void loadSpecificChromsomesOrRegions(
            final CommandLine cmd, final List<String> chromosomes, final List<ChrBaseRegion> regions, final Logger logger) throws ParseException
    {
        if(cmd.hasOption(SPECIFIC_REGIONS))
        {
            regions.addAll(ChrBaseRegion.loadSpecificRegions(cmd));

            for(ChrBaseRegion region : regions)
            {
                logger.info("filtering for specific region: {}", region);

                if(!chromosomes.contains(region.Chromosome))
                    chromosomes.add(region.Chromosome);
            }
        }
        else if(cmd.hasOption(SPECIFIC_CHROMOSOMES))
        {
            logger.info("filtering for specific chromosomes: {}", cmd.getOptionValue(SPECIFIC_CHROMOSOMES));
            chromosomes.addAll(loadSpecificChromsomes(cmd));
        }
    }
}

