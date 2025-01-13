package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class UnmappedRegions
{
    public static final String UNMAP_REGIONS_FILE = "unmap_regions";
    public static final String UNMAP_REGIONS_FILE_DESC = "Genome regions to unmap reads";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(UNMAP_REGIONS_FILE, false, UNMAP_REGIONS_FILE_DESC);
    }

    private static final Logger LOGGER = LogManager.getLogger(UnmappedRegions.class);

    public static Map<String,List<UnmappingRegion>> loadUnmapRegions(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            final String delim = FileDelimiters.inferFileDelimiter(filename);

            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);
            int chrIndex = getChromosomeFieldIndex(fieldIndexMap);
            int posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
            int posEndIndex = getPositionEndFieldIndex(fieldIndexMap);
            int depthIndex = fieldIndexMap.get("MaxDepth");

            Map<String, List<UnmappingRegion>> chrLocationsMap = Maps.newHashMap();

            for(String line : lines)
            {
                String[] values = line.split(delim, -1);

                String chromosome = values[chrIndex];

                List<UnmappingRegion> regions = chrLocationsMap.get(chromosome);
                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    chrLocationsMap.put(chromosome, regions);
                }

                UnmappingRegion lastRegion = (regions.isEmpty()) ? null : regions.get(regions.size() - 1);

                int posStart = Integer.parseInt(values[posStartIndex]);
                int posEnd = Integer.parseInt(values[posEndIndex]);
                int maxDepth = Integer.parseInt(values[depthIndex]);

                UnmappingRegion region = new UnmappingRegion(posStart, posEnd, maxDepth);

                if(lastRegion != null)
                {
                    // too difficult to merge
                    if(lastRegion.overlaps(region))
                    {
                        LOGGER.error("unmap regions overlap: current({}) next({})", lastRegion, region);
                        return null;
                    }
                    else if(region.end() < lastRegion.start())
                    {
                        LOGGER.error("unmap regions not sorted: current({}) next({})", lastRegion, region);
                        return null;
                    }
                }

                regions.add(region);
            }

            return chrLocationsMap;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to read high-depth regions file {}: {}", filename, e.toString());
            System.exit(1);
            return null;
        }
    }
}
