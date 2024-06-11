package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ToFastqUtils
{
    public static final String R1 = "R1";
    public static final String R2 = "R2";
    public static final String UNPAIRED = "Unpaired";

    public static boolean canIgnoreRead(SAMRecord read)
    {
        return read.isSecondaryOrSupplementary() || read.hasAttribute(CONSENSUS_READ_ATTRIBUTE);
    }

    public static SamReader openSamReader(final ToFastqConfig config)
    {
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
        if(config.RefGenomeFile != null)
        {
            samReaderFactory = samReaderFactory.referenceSequence(new File(config.RefGenomeFile));
        }
        return samReaderFactory.open(new File(config.BamFile));
    }

    public static List<SAMReadGroupRecord> getReadGroups(final ToFastqConfig config)
    {
        try(SamReader samReader = ToFastqUtils.openSamReader(config))
        {
            return samReader.getFileHeader().getReadGroups();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    public static void mergeFastqs(Stream<String> inFastqPaths, String outFastq, boolean deleteIfEmpty)
    {
        // whether it is .fastq or .fastq.gz file, we can just simply concat them
        Path outFastqPath = Path.of(outFastq);
        try(OutputStream out = Files.newOutputStream(outFastqPath))
        {
            Iterable<String> inFastqPathsIterable = inFastqPaths::iterator;
            for(String inFastqs : inFastqPathsIterable)
            {
                Path inFastqPath = Path.of(inFastqs);
                if(Files.exists(inFastqPath))
                {
                    Files.copy(inFastqPath, out);
                    Files.delete(inFastqPath);
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        if(deleteIfEmpty)
        {
            try
            {
                if(Files.size(outFastqPath) == 0)
                {
                    Files.delete(outFastqPath);
                }
            }
            catch(IOException e)
            {
                throw new UncheckedIOException(e);
            }
        }
    }

    public static String formFilename(final String filePrefix, final String readType)
    {
        return String.format("%s.%s.fastq.gz", filePrefix, readType);
    }
}
