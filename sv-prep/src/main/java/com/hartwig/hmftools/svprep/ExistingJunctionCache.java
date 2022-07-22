package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.JunctionData;

public class ExistingJunctionCache
{
    private final Map<String,List<JunctionData>> mChrJunctions;

    public ExistingJunctionCache()
    {
        mChrJunctions = Maps.newHashMap();
    }

    public List<JunctionData> getRegionJunctions(final ChrBaseRegion region)
    {
        List<JunctionData> regionJunctions = Lists.newArrayList();

        List<JunctionData> chrJunctions = mChrJunctions.get(region.Chromosome);

        if(chrJunctions == null)
            return regionJunctions;

        // extract and remove matched junctions, assumes they are sequential in the input file
        int index = 0;
        boolean regionFound = false;
        while(index < chrJunctions.size())
        {
            JunctionData junctionData = chrJunctions.get(index);

            if(region.containsPosition(junctionData.Position))
            {
                regionFound = true;
                chrJunctions.remove(index);
                regionJunctions.add(junctionData);
                continue;
            }
            else
            {
                if(regionFound)
                    break;
            }

            ++index;
        }

        return regionJunctions;
    }

    public boolean loadJunctions(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return true;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIM);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int orientIndex = fieldsIndexMap.get("Orientation");

            List<JunctionData> junctionDataList = null;
            String currentChromosome = "";

            int junctionCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIM, -1);

                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                byte orientation = Byte.parseByte(values[orientIndex]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    junctionDataList = Lists.newArrayList();
                    mChrJunctions.put(chromosome, junctionDataList);
                }

                junctionDataList.add(new JunctionData(position, orientation, null));
                ++junctionCount;
            }

            SV_LOGGER.info("loaded {} existing junctions from file: {}", junctionCount, filename);
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read existing junctions file({})", filename, exception.toString());
            return false;
        }

        return true;
    }
}
