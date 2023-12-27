package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class Junction implements Comparable<Junction>
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

    /*
    public final int junctionFragments();
    public final int supportFragments();
    public final int discordantFragments();
    public final int lowMapQualityFragments();
    public final int maxMapQuality();
    public final int maxSoftClipLength();
    public final boolean hasPolyAT();
    public final boolean isIndel();
    public final boolean isHotspot();
    public final String softClipBases();
    public final String initialReadId();
    */

    public Junction(final String chromosome, final int position, final byte orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
    }

    // convenience and poossible temporary
    public String chromosome() { return Chromosome; }
    public int position() { return Position; }
    public Direction orientation() { return Orientation == POS_STRAND ? Direction.FORWARDS : Direction.REVERSE; }

    public String toString()
    {
        return format("%s:%d:%d", Chromosome, Position, Orientation);
    }

    public boolean isLocalMatch(final Junction other)
    {
        return Position == other.Position && Orientation == other.Orientation;
    }

    @Override
    public int compareTo(final Junction other)
    {
        if(!Chromosome.equals(other.Chromosome))
        {
            int firstChrRank = HumanChromosome.chromosomeRank(Chromosome);
            int secondChrRank = HumanChromosome.chromosomeRank(other.Chromosome);

            return firstChrRank < secondChrRank ? -1 : 1;
        }

        if(Position == other.Position)
        {
            if(Orientation == other.Orientation)
                return 0;

            return Orientation == POS_STRAND ? -1 : 1;
        }

        return Position < other.Position ? -1 : 1;
    }

    public static Map<String,List<Junction>> loadJunctions(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return null;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            Map<String,List<Junction>> chrJunctionsMap = Maps.newHashMap();

            String line = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int orientIndex = fieldsIndexMap.get(FLD_ORIENTATION);

            List<Junction> junctionDataList = null;
            String currentChromosome = "";

            int junctionCount = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                byte orientation = Byte.parseByte(values[orientIndex]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    junctionDataList = Lists.newArrayList();
                    chrJunctionsMap.put(chromosome, junctionDataList);
                }

                junctionDataList.add(new Junction(chromosome, position, orientation));
                ++junctionCount;
            }

            SV_LOGGER.info("loaded {} existing junctions from file: {}", junctionCount, filename);

            return chrJunctionsMap;
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read existing junctions file({})", filename, exception.toString());
            return null;
        }
    }

    public static void mergeJunctions(final Map<String,List<Junction>> existingMap, final Map<String,List<Junction>> newMap)
    {
        for(Map.Entry<String,List<Junction>> entry : newMap.entrySet())
        {
            String newChromosome = entry.getKey();
            List<Junction> newJunctions = entry.getValue();

            List<Junction> existingJunctions = existingMap.get(newChromosome);

            if(existingJunctions == null)
            {
                existingJunctions = Lists.newArrayList(newJunctions);
                existingMap.put(newChromosome, existingJunctions);
            }
            else
            {
                for(Junction newJunction : newJunctions)
                {
                    if(existingJunctions.stream().noneMatch(x -> x.isLocalMatch(newJunction)))
                        existingJunctions.add(newJunction);
                }
            }

            Collections.sort(existingJunctions);
        }
    }
}
