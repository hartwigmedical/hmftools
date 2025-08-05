package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class BedFileReader
{
    private static final Logger LOGGER = LogManager.getLogger(BedFileReader.class);

    public static boolean loadBedFile(final String filename, final List<ChrBaseRegion> regions)
    {
        try
        {
            regions.addAll(loadBedFile(filename));
            return true;
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load BED file: {}, error: {}", filename, e.toString());
            return false;
        }
    }

    public static List<ChrBaseRegion> loadBedFile(final String filename) throws Exception
    {
        return loadBedFile(filename, true);
    }

    public static List<ChrBaseRegion> loadBedFile(
            final String filename,
            boolean checkSortedMerged) throws Exception
    {
        return loadBedFile(filename, checkSortedMerged, BedLine.factory());
    }

    public static <T extends ChrBaseRegion>  List<T> loadBedFile(
            final String filename,
            boolean checkSortedMerged,
            Function<String, T> factory) throws IOException
    {
        List<String> lines = FileWriterUtils.readLines(filename);
        return loadBedFile(lines, checkSortedMerged, factory);
    }

    public static List<ChrBaseRegion> loadBedFile(final List<String> lines)
    {
        return loadBedFile(lines, true);
    }

    public static List<ChrBaseRegion> loadBedFile(final List<String> lines, boolean checkSortedMerged)
    {
        return loadBedFile(lines, checkSortedMerged, BedLine.factory());
    }

    public static <T extends ChrBaseRegion>  List<T> loadBedFile(
            final List<String> lines,
            boolean checkSortedMerged,
            Function<String, T> factory)
    {
        List<T> regions = Lists.newArrayList();

        for(String line : lines)
        {
            if(line.contains(FLD_CHROMOSOME))
            {
                continue;
            }
            regions.add(factory.apply(line));
        }

        if(checkSortedMerged)
        {
            ChrBaseRegion.checkMergeOverlaps(regions, true);
        }

        return regions;
    }

    public static Map<Chromosome, List<BaseRegion>> loadBedFileChrMap(final String filename)
    {
        return loadBedFileChrMap(filename, false);
    }

    public static Map<Chromosome, List<BaseRegion>> loadBedFileChrMap(final String filename, boolean checkSortedMerged)
    {
        return loadBedFileChrMap(
                filename,
                checkSortedMerged,
                BedLine.factory(),
                chrBaseRegion -> new BaseRegion(chrBaseRegion.start(), chrBaseRegion.end()));
    }

    public static <R extends BaseRegion, T extends ChrBaseRegion> Map<Chromosome, List<R>> loadBedFileChrMap(
            final String filename,
            boolean checkSortedMerged,
            Function<String, T> factory,
            Function<T, R> converter)
    {
        final Map<Chromosome, List<R>> chrRegionsMap = Maps.newHashMap();

        try
        {
            List<String> lines = FileWriterUtils.readLines(filename);

            List<R> chrRegions = null;
            Chromosome currentChr = null;

            for(String line : lines)
            {
                if(line.contains(FLD_CHROMOSOME))
                {
                    continue;
                }
                T chrRegion = factory.apply(line);
                if(chrRegion == null)
                {
                    continue;
                }

                if(currentChr != chrRegion.humanChromosome())
                {
                    currentChr = chrRegion.humanChromosome();
                    chrRegions = Lists.newArrayList();
                    chrRegionsMap.put(currentChr, chrRegions);
                }

                assert chrRegions != null;
                chrRegions.add(converter.apply(chrRegion));
            }

            if(checkSortedMerged)
            {
                chrRegionsMap.values().forEach(BaseRegion::checkMergeOverlaps);
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load BED file({}): {}", filename, e.toString());
            return null;
        }

        return chrRegionsMap;
    }
}
