package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RegionPerfTracker
{
    public static final int BIN_SIZE = 1_000_000;

    private static final String UNMAPPED_CONTIG = "unmapped";
    private static final double NANOS_IN_SECOND = 1_000_000_000.0;

    private final Map<String, RegionPerf> mRegions;

    public RegionPerfTracker()
    {
        mRegions = new HashMap<>();
    }

    public void add(final String contig, final int position, final long elapsedNanos, final int readCount)
    {
        String chromosome = contig != null ? contig : UNMAPPED_CONTIG;
        int binStart = contig != null ? (position / BIN_SIZE) * BIN_SIZE : 0;

        RegionPerf region = mRegions.computeIfAbsent(
                chromosome + ":" + binStart, key -> new RegionPerf(chromosome, binStart));

        region.add(elapsedNanos, readCount);
    }

    public void merge(final RegionPerfTracker other)
    {
        for(Map.Entry<String, RegionPerf> entry : other.mRegions.entrySet())
        {
            RegionPerf source = entry.getValue();
            RegionPerf region = mRegions.computeIfAbsent(
                    entry.getKey(), key -> new RegionPerf(source.Chromosome, source.PosStart));

            region.merge(source);
        }
    }

    public int regionCount() { return mRegions.size(); }

    public void write(final String filename, final double slowRegionSeconds)
    {
        List<RegionPerf> regions = new ArrayList<>(mRegions.values());
        regions.sort(Comparator.comparingLong((RegionPerf region) -> region.TotalNanos).reversed());

        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            writer.write(String.join(
                    TSV_DELIM, "Chromosome", "PosStart", "PosEnd", "Groups", "Reads", "TotalSeconds", "MaxSeconds"));
            writer.newLine();

            for(RegionPerf region : regions)
            {
                writer.write(String.join(
                        TSV_DELIM,
                        region.Chromosome,
                        String.valueOf(region.PosStart),
                        String.valueOf(region.PosStart + BIN_SIZE - 1),
                        String.valueOf(region.Groups),
                        String.valueOf(region.Reads),
                        String.format("%.3f", region.TotalNanos / NANOS_IN_SECOND),
                        String.format("%.3f", region.MaxNanos / NANOS_IN_SECOND)));

                writer.newLine();
            }

            TARS_LOGGER.info("wrote region perf for {} regions to {}", regions.size(), filename);
        }
        catch(IOException e)
        {
            TARS_LOGGER.warn("failed to write region perf {}: {}", filename, e.toString());
        }

        logSlowRegions(regions, slowRegionSeconds);
    }

    private static void logSlowRegions(final List<RegionPerf> regions, final double slowRegionSeconds)
    {
        for(RegionPerf region : regions)
        {
            double maxSeconds = region.MaxNanos / NANOS_IN_SECOND;
            if(maxSeconds < slowRegionSeconds)
            {
                continue;
            }

            TARS_LOGGER.info("slow region({}:{}-{}) groups({}) reads({}) total({}) max({})",
                    region.Chromosome, region.PosStart, region.PosStart + BIN_SIZE - 1, region.Groups, region.Reads,
                    String.format("%.3f", region.TotalNanos / NANOS_IN_SECOND), String.format("%.3f", maxSeconds));
        }
    }

    private static class RegionPerf
    {
        public final String Chromosome;
        public final int PosStart;

        public long Groups;
        public long Reads;
        public long TotalNanos;
        public long MaxNanos;

        RegionPerf(final String chromosome, final int posStart)
        {
            Chromosome = chromosome;
            PosStart = posStart;
        }

        void add(final long elapsedNanos, final int readCount)
        {
            ++Groups;
            Reads += readCount;
            TotalNanos += elapsedNanos;
            MaxNanos = Math.max(MaxNanos, elapsedNanos);
        }

        void merge(final RegionPerf other)
        {
            Groups += other.Groups;
            Reads += other.Reads;
            TotalNanos += other.TotalNanos;
            MaxNanos = Math.max(MaxNanos, other.MaxNanos);
        }
    }
}
