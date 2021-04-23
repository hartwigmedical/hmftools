package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SPLICE_SITE_FILE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.RegionReadData;

public class SpliceSiteCounter
{
    private final Map<Integer,int[]> mSiteCounts;
    private final BufferedWriter mWriter;

    public static final int SPLICE_SITE_TRAVERSED = 0;
    public static final int SPLICE_SITE_SUPPORT = 1;
    public static final int SPLICE_SITE_SKIPPED = 2;

    public SpliceSiteCounter(final BufferedWriter writer)
    {
        mWriter = writer;
        mSiteCounts = Maps.newHashMap();
    }

    public final Map<Integer,int[]> getSiteCounts() { return mSiteCounts; }
    public void clear() { mSiteCounts.clear(); }

    public void registerSpliceSiteSupport(
            final List<int[]> readMappedCoords1, final List<int[]> readMappedCoords2, final List<RegionReadData> allRegions)
    {
        final Set<Integer> traversedSites = Sets.newHashSet();
        final Set<Integer> supportedSites = Sets.newHashSet();

        registerSpliceSiteSupport(readMappedCoords1, allRegions, traversedSites, supportedSites);
        registerSpliceSiteSupport(readMappedCoords2, allRegions, traversedSites, supportedSites);

        traversedSites.forEach(x -> addCount(x, SPLICE_SITE_TRAVERSED));
        supportedSites.forEach(x -> addCount(x, SPLICE_SITE_SUPPORT));
    }

    private void registerSpliceSiteSupport(
            final List<int[]> readMappedCoords, final List<RegionReadData> allRegions,
            final Set<Integer> traversedSites, final Set<Integer> supportedSites)
    {
        // for each read region (ie unique exon) record if the read supports its splice junction on each side, or skips it
        if(readMappedCoords.size() <= 1)
            return;

        // Set<Integer> skippedSites = Sets.newHashSet();

        for(int i = 0; i < readMappedCoords.size() - 1; ++i)
        {
            int[] mappedCoordLower = readMappedCoords.get(i);
            int[] mappedCoordUpper = readMappedCoords.get(i + 1);
            int junctionLower = mappedCoordLower[SE_END];
            int junctionUpper = mappedCoordUpper[SE_START];

            // record if a region start or end is traversed entirely by a mapped coord junction

            for(RegionReadData region : allRegions)
            {
                boolean isTraversed = false;

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    if(junctionLower < region.Region.Positions[se] && junctionUpper > region.Region.Positions[se])
                    {
                        isTraversed = true;
                        traversedSites.add(region.Region.Positions[se]);
                    }
                }

                if(isTraversed)
                    continue;

                // does the read match the pre-region boundary on any transcript
                if(region.getPreRegions().stream().anyMatch(x -> x.Region.end() == junctionLower))
                {
                    if(junctionUpper == region.start())
                        supportedSites.add(region.start());
                    // else if(junctionUpper > region.start())
                    //    skippedSites.add(region.start());
                }
                else if(region.getPostRegions().stream().anyMatch(x -> x.start() == junctionUpper))
                {
                    if(junctionLower == region.end())
                        supportedSites.add(region.end());
                    // else if(junctionLower < region.end())
                    //    skippedSites.add(region.end());
                }
            }
        }
    }

    private void addCount(int position, int type)
    {
        int[] counts = mSiteCounts.get(position);

        if(counts == null)
        {
            counts = new int[SPLICE_SITE_SKIPPED + 1];
            mSiteCounts.put(position, counts);
        }

        ++counts[type];
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        if(config.OutputDir == null || config.SampleId == null)
            return null;

        try
        {
            final String outputFileName = config.formOutputFile(SPLICE_SITE_FILE);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneSetId,Chromosome,SpliceSitePosition,TraverseFrags,SupportFrags");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site file: {}", e.toString());
            return null;
        }
    }

    public void writeSpliceSiteData(final GeneCollection geneCollection)
    {
        if(mWriter == null)
            return;

        writeSpliceSiteData(mWriter, geneCollection, mSiteCounts);
    }

    private synchronized static void writeSpliceSiteData(
            final BufferedWriter writer, final GeneCollection geneCollection, final Map<Integer,int[]> siteCounts)
    {
        if(writer == null)
            return;

        try
        {
            for(Map.Entry<Integer,int[]> entry : siteCounts.entrySet())
            {
                int spliceSite = entry.getKey();
                final int[] counts = entry.getValue();

                writer.write(String.format("%s,%s,%d", geneCollection.chrId(), geneCollection.chromosome(), spliceSite));

                writer.write(String.format(",%d,%d", counts[SPLICE_SITE_TRAVERSED], counts[SPLICE_SITE_SUPPORT]));

                // StringJoiner transIdStr = new StringJoiner(SUB_ITEM_DELIM);
                // region.getTransExonRefs().forEach(x -> transIdStr.add(String.valueOf(x.TransId)));

                // writer.write(String.format(",%s", transIdStr.toString()));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site data file: {}", e.toString());
        }
    }
}
