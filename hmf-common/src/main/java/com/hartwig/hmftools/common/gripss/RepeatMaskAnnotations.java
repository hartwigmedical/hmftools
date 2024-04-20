package com.hartwig.hmftools.common.gripss;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RepeatMaskAnnotations
{
    private final Map<String,List<RepeatMaskData>> mChrDataMap;

    public static final String REPEAT_MASK_FILE = "repeat_mask_file";
    public static final String REPEAT_MASK_FILE_DESC = "Repeat mask definitions file";

    private static final Logger LOGGER = LogManager.getLogger(RepeatMaskAnnotations.class);

    public RepeatMaskAnnotations()
    {
        mChrDataMap = Maps.newHashMap();
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(REPEAT_MASK_FILE, false, REPEAT_MASK_FILE_DESC);
    }

    public boolean hasData() { return !mChrDataMap.isEmpty(); }

    public List<RepeatMaskData> findMatches(final String chromosome, final BaseRegion region)
    {
        return findMatches(chromosome, region.start(), region.end());
    }

    public List<RepeatMaskData> findMatches(final ChrBaseRegion region)
    {
        return findMatches(region.Chromosome, region.start(), region.end());
    }

    private static final int SMALL_LINEAR_SEARCH_SIZE = 50;

    public List<RepeatMaskData> findMatches(final String chromosome, final int regionStart, final int regionEnd)
    {
        List<RepeatMaskData> regions = mChrDataMap.get(chromosome);

        if(regions == null)
            return Lists.newArrayList();

        if(regions.size() < SMALL_LINEAR_SEARCH_SIZE)
        {
            return regions.stream()
                    .filter(x -> positionsOverlap(x.Region.start(), x.Region.end(), regionStart, regionEnd)).collect(Collectors.toList());
        }

        // use a binary search since the number of entries is typically > 100K per chromosome
        int currentIndex = regions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = regions.size() - 1;

        List<RepeatMaskData> matchedRegions = Lists.newArrayList();

        while(true)
        {
            RepeatMaskData rmData = regions.get(currentIndex);

            if(regionEnd < rmData.Region.start())
            {
                if(lowerIndex + 1 == currentIndex)
                    break;

                upperIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else if(regionStart > rmData.Region.end())
            {
                // search higher
                if(currentIndex + 1 == upperIndex)
                    break;

                lowerIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else if(positionsOverlap(rmData.Region.start(), rmData.Region.end(), regionStart, regionEnd))
            {
                matchedRegions.add(rmData);
                break;
            }
        }

        // check up and down for further overlaps
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? currentIndex + 1 : currentIndex - 1;

            while(index >= 0 && index < regions.size())
            {
                RepeatMaskData rmData = regions.get(index);

                if(!positionsOverlap(rmData.Region.start(), rmData.Region.end(), regionStart, regionEnd))
                    break;

                matchedRegions.add(rmData);

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }

        return matchedRegions;
    }

    public boolean load(final String filename, final RefGenomeVersion refGenomeVersion)
    {
        if(filename == null)
            return true;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String line = null;
            String currentChr = "";
            List<RepeatMaskData> entries = null;
            int index = 0;

            // first 3 lines contain the header, then expect columns as:
            // SW     perc perc perc  query      position in query           matching       repeat              position in  repeat
            // score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID
            // 0      1    2    3    4           5       6     7           8  9              10                       11 12     13       14
            // 1504   1.3  0.4  1.3  chr1        10001   10468 (249240153) +  (CCCTAA)n      Simple_repeat            1  463    (0)      1
            fileReader.readLine();
            fileReader.readLine();
            fileReader.readLine();

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.trim().split("\\s{1,}", -1);

                String chromosome = refGenomeVersion.versionedChromosome(values[4]);

                if(!HumanChromosome.contains(chromosome))
                    continue;

                if(!chromosome.equals(currentChr))
                {
                    currentChr = chromosome;
                    entries = Lists.newArrayList();
                    mChrDataMap.put(chromosome, entries);
                }

                try
                {
                    // note BED start position adjustment
                    BaseRegion region = new BaseRegion(Integer.parseInt(values[5]) + 1, Integer.parseInt(values[6]));
                    int id = Integer.parseInt(values[14]);
                    int swScore = Integer.parseInt(values[0]);
                    char orientation = values[8].charAt(0);
                    String classType = values[10];
                    String repeat = values[9];

                    entries.add(new RepeatMaskData(id, region, swScore, orientation, repeat, classType));
                }
                catch(Exception e)
                {
                    LOGGER.error("invalid RM file entry: index({}) line({})", line, index);
                    return false;
                }

                ++index;
            }

            LOGGER.info("loaded {} repeat-mask entries from file({})",
                    mChrDataMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load repeat-mask data from file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }
}
