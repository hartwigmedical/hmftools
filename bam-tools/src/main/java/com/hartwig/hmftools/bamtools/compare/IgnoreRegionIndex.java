package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

// per-chromosome sorted index of ignore regions, queried by alignment overlap. Multiple TSV files are
// merged into a single index. Each TSV must have a header row including Chromosome / PosStart / PosEnd
// (extra columns are tolerated). An alignment "hits" the index when any single ignore region overlaps
// more than half the alignment's reference span.
public class IgnoreRegionIndex
{
    private final Map<String, List<BaseRegion>> mRegionsByChromosome;
    private final int mRegionCount;

    private IgnoreRegionIndex(final Map<String, List<BaseRegion>> regionsByChromosome, final int regionCount)
    {
        mRegionsByChromosome = regionsByChromosome;
        mRegionCount = regionCount;
    }

    public int regionCount() { return mRegionCount; }

    public int chromosomeCount() { return mRegionsByChromosome.size(); }

    // returns true if any ignore region on the alignment's chromosome overlaps it by more than half its
    // reference span. Span is (alignmentEnd - alignmentStart + 1); for an unmapped record callers should
    // skip this check entirely (no coords).
    public boolean overlapsAboveHalf(final String chromosome, final int alignmentStart, final int alignmentEnd)
    {
        final List<BaseRegion> regions = mRegionsByChromosome.get(chromosome);
        if(regions == null || regions.isEmpty())
            return false;

        final int span = alignmentEnd - alignmentStart + 1;
        if(span <= 0)
            return false;

        final int threshold = span / 2;

        for(final BaseRegion region : regions)
        {
            if(region.start() > alignmentEnd)
                break;
            if(region.end() < alignmentStart)
                continue;

            final int overlapStart = Math.max(region.start(), alignmentStart);
            final int overlapEnd = Math.min(region.end(), alignmentEnd);
            final int overlap = overlapEnd - overlapStart + 1;

            if(overlap > threshold)
                return true;
        }

        return false;
    }

    public static IgnoreRegionIndex load(final List<String> filenames)
    {
        final Map<String, List<BaseRegion>> merged = new HashMap<>();
        int totalRegions = 0;

        for(final String filename : filenames)
        {
            totalRegions += loadOne(filename, merged);
        }

        for(final List<BaseRegion> regions : merged.values())
        {
            regions.sort(Comparator.comparingInt(BaseRegion::start));
        }

        BT_LOGGER.info("ignore-region index: loaded {} regions across {} chromosomes from {} file(s)",
                totalRegions, merged.size(), filenames.size());

        return new IgnoreRegionIndex(merged, totalRegions);
    }

    private static int loadOne(final String filename, final Map<String, List<BaseRegion>> target)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            if(lines.isEmpty())
                return 0;

            final String delim = FileDelimiters.inferFileDelimiter(filename);
            final String header = lines.get(0);

            final Map<String, Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);
            final int chrIndex = getChromosomeFieldIndex(fieldIndexMap);
            final int posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
            final int posEndIndex = getPositionEndFieldIndex(fieldIndexMap);

            int added = 0;
            for(int i = 1; i < lines.size(); ++i)
            {
                final String line = lines.get(i);
                if(line.isEmpty())
                    continue;

                final String[] values = line.split(delim, -1);
                final String chromosome = values[chrIndex];
                final int posStart = Integer.parseInt(values[posStartIndex]);
                final int posEnd = Integer.parseInt(values[posEndIndex]);

                target.computeIfAbsent(chromosome, k -> new java.util.ArrayList<>()).add(new BaseRegion(posStart, posEnd));
                ++added;
            }

            BT_LOGGER.debug("ignore-region file {}: {} regions", filename, added);
            return added;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read ignore-region file: " + filename, e);
        }
    }
}
