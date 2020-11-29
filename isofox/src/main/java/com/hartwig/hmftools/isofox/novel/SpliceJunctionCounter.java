package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionReadData.SPLICE_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.SPLICE_JUNCTION_TOTAL;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SPLICE_SITE_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUB_ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class SpliceJunctionCounter
{
    public SpliceJunctionCounter()
    {

    }

    public static void assignSpliceJunctionSupport(final List<int[]> readMappedCoords, final List<RegionReadData> regions)
    {
        // for each read region (ie unique exon) record if the read supports its splice junction on each side, or skips it
        if(readMappedCoords.size() <= 1)
            return;

        for(int i = 0; i < readMappedCoords.size() - 1; ++i)
        {
            int[] mappedCoordLower = readMappedCoords.get(i);
            int[] mappedCoordUpper = readMappedCoords.get(i + 1);
            int junctionLower = mappedCoordLower[SE_END];
            int junctionUpper = mappedCoordUpper[SE_START];

            for(RegionReadData region : regions)
            {
                // does the read match the pre-region boundary on any transcript
                if(region.getPreRegions().stream().anyMatch(x -> x.Region.end() == junctionLower))
                {
                    boolean spliceSiteSupported = region.start() == junctionUpper;
                    region.addSpliceJunctionSupport(SE_START, spliceSiteSupported);
                }
                else if(region.getPostRegions().stream().anyMatch(x -> x.start() == junctionUpper))
                {
                    boolean spliceSiteSupported = region.end() == junctionLower;
                    region.addSpliceJunctionSupport(SE_END, spliceSiteSupported);
                }
            }
        }
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        if(config.OutputDir == null || config.SampleId == null)
            return null;

        try
        {
            final String outputFileName = config.formOutputFile(SPLICE_SITE_FILE);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneSetId,Chromosome,RegionStart,RegionEnd");
            writer.write(",StartSupportFrags,StartTotalFrags,EndSupportFrags,EndTotalFrags,TransIds");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeSpliceSiteData(
            final BufferedWriter writer, final GeneCollection geneCollection, final RegionReadData region)
    {
        if(writer == null)
            return;

        final int[][] spliceCounts = region.getSpliceJunctionSupport();

        // skip if this junction has no fragment support of any kind
        if(spliceCounts[SE_START][SPLICE_JUNCTION_TOTAL] == 0 && spliceCounts[SE_END][SPLICE_JUNCTION_TOTAL] == 0)
            return;

        try
        {
            writer.write(String.format("%s,%d,%d", geneCollection.chrId(), region.start(), region.end()));


            writer.write(String.format(",%d,%d,%d,%d",
                    spliceCounts[SE_START][SPLICE_JUNCTION_SUPPORT], spliceCounts[SE_START][SPLICE_JUNCTION_TOTAL],
                    spliceCounts[SE_END][SPLICE_JUNCTION_SUPPORT], spliceCounts[SE_END][SPLICE_JUNCTION_TOTAL]));

            StringJoiner transIdStr = new StringJoiner(SUB_ITEM_DELIM);
            region.getTransExonRefs().forEach(x -> transIdStr.add(String.valueOf(x.TransId)));

            writer.write(String.format(",%s", transIdStr.toString()));
            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site data file: {}", e.toString());
        }
    }
}
