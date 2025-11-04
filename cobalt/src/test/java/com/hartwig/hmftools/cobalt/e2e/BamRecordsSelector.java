package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamRecordsSelector
{
    public static void main(final String[] args)
    {
        File inputBam = new File("/Users/timlavers/work/data/COLO829v003/COLO829v003T.lane_merged.v38.bam");
        File outputDirectory = new File("/Users/timlavers/work/junk");
        File outputBam = new File(outputDirectory, "selected.bam");
        int start = 34491000;
        int end = 34493000;
        new BamRecordsSelector(inputBam, outputBam, "chr2", start, end).run();
    }

    private final File InputBam;
    private final File OutputBam;
    final private String Contig;
    final private int Start;
    final private int End;

    public BamRecordsSelector(final File inputBam, final File outputFile, final String contig, final int start, final int end)
    {
        InputBam = inputBam;
        OutputBam = outputFile;
        Contig = contig;
        Start = start;
        End = end;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();
        CB_LOGGER.info("starting BAM records selector");
        ChrBaseRegion region = new ChrBaseRegion(Contig, Start, End);
        try(SamReader samReader = SamReaderFactory.makeDefault().open(InputBam))
        {
            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();
            try(SAMFileWriter bamWriter = new FastBamWriter(newHeader, OutputBam.getAbsolutePath()))
            {
                BamSlicer bamSlicer = new BamSlicer(0, true, false, false);
                bamSlicer.slice(samReader, region, bamWriter::addAlignment);
            }
            catch(Exception e)
            {
                CB_LOGGER.error("Error processing or writing records", e);
                throw new RuntimeException(e);
            }

            CB_LOGGER.debug("BAM writer closed");
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        CB_LOGGER.info("Selection complete, mins({})", runTimeMinsStr(startTimeMs));
    }
}
