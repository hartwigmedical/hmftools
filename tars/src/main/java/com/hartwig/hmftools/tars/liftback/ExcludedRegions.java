package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Genomic regions reads are excluded from (e.g. RNA rRNA / 7SL / acrocentric contamination zones). Loaded from a
// regions TSV (Chromosome/PosStart/PosEnd) and queried post-lift on lifted genomic coords (excludes()): a primary
// lifting in is unmapped, a supp dropped. Post-lift because a tx-contig read's input coords are chrN_tx and can't
// be tested against the genomic region list; some excluded zones (acrocentric arms) ARE in the transcriptome.
// Per-chromosome regions are sorted for an O(log N) overlap test.
public class ExcludedRegions
{
    private final Map<String, List<ChrBaseRegion>> mRegionsByChr;

    ExcludedRegions(final Map<String, List<ChrBaseRegion>> regionsByChr)
    {
        mRegionsByChr = regionsByChr;
        for(final List<ChrBaseRegion> regions : mRegionsByChr.values())
            regions.sort(Comparator.comparingInt(ChrBaseRegion::start));
    }

    public static ExcludedRegions load(final String filename)
    {
        return new ExcludedRegions(loadRegions(filename));
    }

    // TSV with Chromosome/PosStart/PosEnd columns -> per-chromosome region lists.
    private static Map<String, List<ChrBaseRegion>> loadRegions(final String filename)
    {
        final Map<String, List<ChrBaseRegion>> regionsByChr = Maps.newHashMap();

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            final String delim = inferFileDelimiter(filename);
            final Map<String, Integer> fieldIndexMap = createFieldsIndexMap(lines.get(0), delim);
            final int chrIndex = getChromosomeFieldIndex(fieldIndexMap);
            final int posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
            final int posEndIndex = getPositionEndFieldIndex(fieldIndexMap);

            for(int i = 1; i < lines.size(); ++i)
            {
                final String[] values = lines.get(i).split(delim, -1);
                final String chromosome = values[chrIndex];
                final int posStart = Integer.parseInt(values[posStartIndex]);
                final int posEnd = Integer.parseInt(values[posEndIndex]);
                regionsByChr.computeIfAbsent(chromosome, x -> Lists.newArrayList())
                        .add(new ChrBaseRegion(chromosome, posStart, posEnd));
            }
        }
        catch(IOException e)
        {
            TARS_LOGGER.error("failed to read excluded regions file {}: {}", filename, e.toString());
            return null;
        }

        return regionsByChr;
    }

    // true if [posStart, posEnd] on chromosome overlaps an excluded region. Called post-lift on genomic coords.
    public boolean excludes(final String chromosome, final int posStart, final int posEnd)
    {
        return overlaps(chromosome, posStart, posEnd);
    }

    private boolean overlaps(final String chromosome, final int posStart, final int posEnd)
    {
        final List<ChrBaseRegion> regions = mRegionsByChr.get(chromosome);
        if(regions == null || regions.isEmpty())
            return false;

        // rightmost region whose start <= posEnd; non-overlapping regions, so it's the only candidate.
        int lo = 0;
        int hi = regions.size() - 1;
        int candidate = -1;
        while(lo <= hi)
        {
            final int mid = (lo + hi) >>> 1;
            if(regions.get(mid).start() <= posEnd)
            {
                candidate = mid;
                lo = mid + 1;
            }
            else
            {
                hi = mid - 1;
            }
        }

        return candidate >= 0 && regions.get(candidate).end() >= posStart;
    }
}
