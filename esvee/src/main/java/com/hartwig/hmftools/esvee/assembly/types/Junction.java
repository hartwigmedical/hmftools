package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.fromByteStr;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.compareJunctions;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_EXACT_SUPPORT_FRAGS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_EXTRA_INFO;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_HOTSPOT_JUNCTION;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_INDEL_JUNCTION;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_JUNCTION_FRAGS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_OTHER_SUPPORT_FRAGS;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public class Junction implements Comparable<Junction>
{
    public final String Chromosome;
    public final int Position;
    public final Orientation Orient;

    public final boolean DiscordantOnly;
    public final boolean Hotspot;

    private boolean mIndelBased;
    public final String mDetails;
    private int mRawDiscordantPosition;
    private IndelCoords mIndelCoords; // the consensus indel if applicable

    public Junction(final String chromosome, final int position, final Orientation orientation)
    {
        this(chromosome, position, orientation, false, false, false);
        mRawDiscordantPosition = -1;
    }

    public Junction(
            final String chromosome, final int position, final Orientation orientation, final boolean discordantOnly,
            final boolean indelBased, final boolean hotspot)
    {
        Chromosome = chromosome;
        Position = position;
        Orient = orientation;
        DiscordantOnly = discordantOnly;
        mIndelBased = indelBased;
        Hotspot = hotspot;

        if(DiscordantOnly || mIndelBased || Hotspot)
        {
            StringJoiner sj = new StringJoiner("/");
            if(DiscordantOnly)
                sj.add("disc-only");
            if(mIndelBased)
                sj.add("indel");
            if(Hotspot)
                sj.add("hotspot");

            mDetails = sj.toString();
        }
        else
        {
            mDetails = "";
        }

        mIndelCoords = null;
    }

    public boolean isForward() { return Orient.isForward(); }
    public boolean isReverse() { return Orient.isReverse(); }

    public void markAsIndel() { mIndelBased = true; }
    public boolean indelBased() { return mIndelBased; }

    public void setIndelCoords(final IndelCoords indelCoords) { mIndelCoords = indelCoords; }
    public IndelCoords indelCoords() { return mIndelCoords; }

    public int rawDiscordantPosition() { return mRawDiscordantPosition; }
    public void setRawDiscordantPosition(int position) { mRawDiscordantPosition = position; }

    public String toString()
    {
        if(DiscordantOnly || mIndelBased || Hotspot)
        {
            return format("%s:%d:%d %s",Chromosome, Position, Orient.asByte(), mDetails);
        }

        return format("%s:%d:%d", Chromosome, Position, Orient.asByte());
    }

    // for display and logging
    public String coords() { return format("%s:%d:%d", Chromosome, Position, Orient.asByte()); }
    public String coordsTyped() { return coordsTyped(false); }

    public String coordsTyped(boolean useRawDiscordant)
    {
        if(DiscordantOnly)
            return format("%s:%d:%d:D", Chromosome, useRawDiscordant ? mRawDiscordantPosition : Position, Orient.asByte());
        else if(indelBased())
            return format("%s:I", coords());
        else
            return coords();
    }

    public boolean isLocalMatch(final Junction other)
    {
        return Position == other.Position && Orient == other.Orient;
    }

    public boolean lower(final Junction other)
    {
        return Position < other.Position || (Position == other.Position && Orient.isForward() && other.Orient.isReverse());
    }

    public boolean higher(final Junction other) { return !lower(other); }

    @Override
    public int compareTo(final Junction other)
    {
        return compareJunctions(Chromosome, other.Chromosome, Position, other.Position, Orient, other.Orient);
    }

    public static Map<String,List<Junction>> loadJunctions(
            final String filename, final SpecificRegions specificRegions, int minJunctionFrags, int minHotspotFrags, int minDiscordantFrags)
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

            Integer juncFragsIndex = fieldsIndexMap.get(FLD_JUNCTION_FRAGS);
            Integer otherJuncFragsIndex = fieldsIndexMap.get(FLD_EXACT_SUPPORT_FRAGS);

            Integer otherSupportFragsIndex = fieldsIndexMap.containsKey(FLD_OTHER_SUPPORT_FRAGS) ?
                    fieldsIndexMap.get(FLD_OTHER_SUPPORT_FRAGS) : fieldsIndexMap.get("DiscordantFrags"); // old name

            Integer indelIndex = fieldsIndexMap.get(FLD_INDEL_JUNCTION);
            Integer hotspotIndex = fieldsIndexMap.get(FLD_HOTSPOT_JUNCTION);
            Integer extraInfoIndex = fieldsIndexMap.get(FLD_EXTRA_INFO);

            List<Junction> junctionDataList = null;
            String currentChromosome = "";

            int junctionCount = 0;
            int rawJunctionCount = 0;
            int discordantCount = 0;
            int indelCount = 0;
            int hotspotCount = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                Orientation orientation = fromByteStr(values[orientIndex]);

                if(!specificRegions.includeChromosome(chromosome))
                    continue;

                if(!specificRegions.includePosition(chromosome, position))
                    continue;

                ++rawJunctionCount;

                int junctionFrags = juncFragsIndex != null ? Integer.parseInt(values[juncFragsIndex]) : 0;
                int otherJunctionFrags = otherSupportFragsIndex != null ? Integer.parseInt(values[otherJuncFragsIndex]) : 0;
                int otherSupportFrags = otherSupportFragsIndex != null ? Integer.parseInt(values[otherSupportFragsIndex]) : 0;

                boolean discordantOnly = junctionFrags == 0 && otherSupportFrags > 0;
                boolean indel = indelIndex != null && Boolean.parseBoolean(values[indelIndex]);
                boolean hotspot = hotspotIndex != null && Boolean.parseBoolean(values[hotspotIndex]);

                if(hotspot)
                {
                    if(junctionFrags < minHotspotFrags)
                        continue;

                    ++hotspotCount;
                }
                else if(discordantOnly)
                {
                    int maxRemoteFrags = extraInfoIndex != null ? Integer.parseInt(values[extraInfoIndex]) : otherSupportFrags;

                    if(maxRemoteFrags < minDiscordantFrags)
                        continue;

                    ++discordantCount;
                }
                else
                {
                    if(junctionFrags + otherJunctionFrags < minJunctionFrags)
                        continue;

                    if(indel)
                        ++indelCount;
                }

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    junctionDataList = Lists.newArrayList();
                    chrJunctionsMap.put(chromosome, junctionDataList);
                }

                junctionDataList.add(new Junction(chromosome, position, orientation, discordantOnly, indel, hotspot));
                ++junctionCount;
            }

            int splitCount = junctionCount - discordantCount - indelCount;
            SV_LOGGER.info("loaded {} junctions (raw={}) types(split={} indel={} discordant={} hotspot={}) from file: {}",
                    junctionCount, rawJunctionCount, splitCount, indelCount, discordantCount, hotspotCount, filename);

            chrJunctionsMap.values().forEach(x -> Collections.sort(x));

            return chrJunctionsMap;
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read junctions file({})", filename, exception.toString());
            return null;
        }
    }

    public static Junction fromConfigStr(final String configStr)
    {
        String[] values = configStr.split(":", -1);
        if(values.length < 3)
            return null;

        String chromosome = values[0];
        int position = Integer.parseInt(values[1]);
        Orientation orientation = fromByteStr(values[2]);

        boolean discordantOnly = false;
        boolean indelBased = false;

        if(values.length >= 4)
        {
            if(values[3].equals("D"))
                discordantOnly = true;
            else if(values[3].equals("I"))
                indelBased = true;
        }

        return new Junction(chromosome, position, orientation, discordantOnly, indelBased, false);
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
                        moveExisting = true;
                    }
                    else if(newJunction.higher(existingJunction))
                    {
                        moveExisting = true;
                    }
                    else if(newJunction.lower(existingJunction))
                    {
                        existingJunctions.add(existingIndex, newJunction);
                        ++existingIndex;
                        moveNew = true;
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

    public static boolean validateJunctionMap(final Map<String,List<Junction>> map)
    {
        // checks order and for duplicates
        for(List<Junction> junctions : map.values())
        {
            for(int i = 0; i < junctions.size() - 1; ++i)
            {
                Junction junction = junctions.get(i);
                Junction nextJunction = junctions.get(i + 1);

                if(junction.isLocalMatch(nextJunction))
                {
                    SV_LOGGER.error("junction({}:{}) duplicate", i, junction);
                    return false;
                }

                if(junction.Position > nextJunction.Position)
                {
                    SV_LOGGER.error("junction({}:{}) after next({})", i, junction, nextJunction);
                    return false;
                }

                // duplicates
                for(int j = i + 1; j < junctions.size(); ++j)
                {
                    if(junction.isLocalMatch(junctions.get(j)))
                    {
                        SV_LOGGER.error("junction({}:{}) duplicate with later index({})", i, junction, j);
                        return false;
                    }
                }
            }
        }

        return true;
    }
}
