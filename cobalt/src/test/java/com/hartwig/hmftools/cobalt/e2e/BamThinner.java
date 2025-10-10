package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Random;

import com.hartwig.hmftools.common.bam.FastBamWriter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamThinner
{
    private final File InputBam;
    private final File OutputDirectory;
    private final RecordMatcher Matcher;
    private final Random mRandom = new Random();
    private final double FractionOfRecordsToKeep;

    public BamThinner(final File inputBam, final File outputDirectory, final double fractionOfRecordsToKeep)
    {
        InputBam = inputBam;
        OutputDirectory = outputDirectory;
        Matcher = new RecordMatcher(fractionOfRecordsToKeep);
        FractionOfRecordsToKeep = fractionOfRecordsToKeep;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();
        CB_LOGGER.info("starting BAM thinner");

        try(SamReader samReader = SamReaderFactory.makeDefault().open(InputBam))
        {
            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();

            File interimOutputFile = new File(OutputDirectory, "thinned_" + InputBam.getName());
            try(SAMFileWriter bamWriter = new FastBamWriter(newHeader, interimOutputFile.getAbsolutePath()))
            {
//                samReader.forEach(record -> Matcher.process(record).forEach(bamWriter::addAlignment));
                samReader.forEach(samRecord -> {
                    if (mRandom.nextDouble() <= FractionOfRecordsToKeep) {
                        bamWriter.addAlignment(samRecord);
                    }
                });
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
        CB_LOGGER.info("Thinning complete, mins({})", runTimeMinsStr(startTimeMs));
        CB_LOGGER.info("Number of records seen ({})", Matcher.processed());
        CB_LOGGER.info("Number of records retained ({})", Matcher.retained());
    }

    public static void main(final String[] args)
    {
        File inputBam = new File("/Users/timlavers/work/data/COLO829/COLO829R.bam");
//        File inputBam = new File("/Users/timlavers/work/scratch/datasets/pmhaem/bam/Sample_13927535.bam");
        File outputDirectory = new File("/Users/timlavers/work/junk/rubbish");
        new BamThinner(inputBam, outputDirectory, 0.1).run();
    }
}
