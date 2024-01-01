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
    public Direction direction() { return Orientation == POS_STRAND ? Direction.FORWARDS : Direction.REVERSE; }

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

            SV_LOGGER.info("loaded {} junctions from file: {}", junctionCount, filename);

            return chrJunctionsMap;
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read junctions file({})", filename, exception.toString());
            return null;
        }
    }

    public static void mergeJunctions(final Map<String,List<Junction>> existingMap, final Map<String,List<Junction>> newMap)
    {
        if(existingMap.isEmpty())
        {
            existingMap.putAll(newMap);
            return;
        }

        for(Map.Entry<String,List<Junction>> entry : newMap.entrySet())
        {
            List<Junction> newJunctions = entry.getValue();

            if(newJunctions.isEmpty())
                continue;

            String newChromosome = entry.getKey();

            List<Junction> existingJunctions = existingMap.get(newChromosome);

            if(existingJunctions == null || existingJunctions.isEmpty())
            {
                existingMap.put(newChromosome, Lists.newArrayList(newJunctions));
                continue;
            }

            // walk through together
            int existingIndex = 0;
            int newIndex = 0;

            Junction existingJunction = existingJunctions.get(existingIndex);
            Junction newJunction = newJunctions.get(newIndex);

            while(true)
            {
                boolean moveNew = false;
                boolean moveExisting = false;

                if(existingJunction == null)
                {
                    existingJunctions.add(newJunction);
                    moveNew = true;
                }
                else
                {
                    if(existingJunction.isLocalMatch(newJunction))
                    {
                        // a match so move both, and without inserting a new entry
                        moveNew = true;

                        // preempt a match with the next but out of order in orientation
                        if(newIndex < newJunctions.size() - 1 && existingIndex > 0
                        && existingJunctions.get(existingIndex - 1).isLocalMatch(newJunctions.get(newIndex + 1)))
                        {
                            --existingIndex;
                            existingJunction = existingJunctions.get(existingIndex);
                        }
                        else
                        {
                            moveExisting = true;
                        }
                    }
                    else if(newJunction.Position > existingJunction.Position)
                    {
                        moveExisting = true;
                    }
                    else if(newJunction.Position < existingJunction.Position)
                    {
                        existingJunctions.add(existingIndex, newJunction);
                        ++existingIndex;
                        moveNew = true;
                    }
                    else if(existingJunction.Orientation != newJunction.Orientation)
                    {
                        moveExisting = true;
                    }
                    else
                    {
                        SV_LOGGER.error("junction merge failed: existing({}:{}) vs new({}:{})",
                                existingIndex, existingJunction, newIndex, newJunction);
                        return;
                    }
                }

                if(moveNew)
                {
                    ++newIndex;

                    if(newIndex >= newJunctions.size())
                        break;

                    newJunction = newJunctions.get(newIndex);
                }

                if(moveExisting)
                {
                    ++existingIndex;

                    if(existingIndex >= existingJunctions.size())
                        existingJunction = null;
                    else
                        existingJunction = existingJunctions.get(existingIndex);
                }
            }
        }
    }
}
