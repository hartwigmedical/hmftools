package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.*;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

public class AltContigRemapper
{
    private final AltContigRemapperConfig mConfig;

    public AltContigRemapper(final AltContigRemapperConfig config)
    {
        mConfig = config;
    }

    static boolean isRelevantDictionaryItem(SAMSequenceRecord samSequenceRecord)
    {
        return !samSequenceRecord.getSequenceName().toLowerCase().startsWith("hla");
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();
        BT_LOGGER.info("starting alt contig remapper");
        try(SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.mOrigBamFile)))
        {

            // The header in the rewritten file needs to be the same
            // as the header in the original file but with the hla dictionary items removed.
            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();
            SAMSequenceDictionary dictionaryWithHlaAltsRemoved = new SAMSequenceDictionary();
            fileHeader.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(AltContigRemapper::isRelevantDictionaryItem)
                    .forEach(dictionaryWithHlaAltsRemoved::addSequence);
            newHeader.setSequenceDictionary(dictionaryWithHlaAltsRemoved);

            // The initial output is unsorted.
            String interimOutputFileName = mConfig.mOutputFile + ".unsorted";
            File interimOutputFile = new File(interimOutputFileName);
            try(SAMFileWriter bamWriter = new FastBamWriter(newHeader, interimOutputFileName))
            {
                BT_LOGGER.info("New BAM writer created.");

                final BwaHlaRecordAligner aligner = new BwaHlaRecordAligner(mConfig.aligner(), newHeader, mConfig.mRefGenVersion);
                HlaTransformer transformer = new HlaTransformer(aligner);
                samReader.forEach(record -> transformer.process(record).forEach(bamWriter::addAlignment));

                BT_LOGGER.info("Input file processed. Number of HLA records: " + transformer.numberOfHlaRecordsProcessed());
                // Deal with any unmatched reads.
                // Don't map these - log an error and write them out as they are
                List<SAMRecord> unmatched = transformer.unmatchedRecords();
                if(!unmatched.isEmpty())
                {
                    BT_LOGGER.warn("Some HLA contig records were unmatched. " + unmatched);
                    unmatched.forEach(bamWriter::addAlignment);
                }
                else
                {
                    BT_LOGGER.info("No HLA contig records were unmatched.");
                }
            }
            catch(Exception e)
            {
                BT_LOGGER.error("Error processing or writing records.", e);
                throw new RuntimeException(e);
            }
            BT_LOGGER.info("BAM Writer closed.");

            // If the samtools path has been provided, sort the output. Else simply rename the unsorted file.
            if(mConfig.mBamToolPath != null)
            {
                BT_LOGGER.info("Sorting output. Threads: " + mConfig.mThreads);
                writeSortedBam(interimOutputFileName, mConfig.mOutputFile, mConfig.mBamToolPath, mConfig.mThreads);
                BT_LOGGER.info("Sorting complete.");
            }
            else
            {
                File outputFile = new File(mConfig.mOutputFile);
                boolean renamed = interimOutputFile.renameTo(outputFile);
                if(!renamed)
                {
                    BT_LOGGER.warn("Could not rename " + interimOutputFile + " to " + outputFile);
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        BT_LOGGER.info("Remapping complete in {} minutes.", runTimeMinsStr(startTimeMs));
    }

    private static void writeSortedBam(final String unsortedBam, final String sortedBam, final String bamToolPath, final int threads)
    {
        if(bamToolPath == null)
        {
            return;
        }

        BT_LOGGER.info("writing sorted BAM: {}", sortedBam);

        BamToolName toolName = fromPath(bamToolPath);

        boolean success = BamOperations.sortBam(toolName, bamToolPath, unsortedBam, sortedBam, threads);

        if(success && toolName == BamToolName.SAMTOOLS)
        {
            success = BamOperations.indexBam(toolName, bamToolPath, sortedBam, threads);
        }

        if(success)
        {
            try
            {
                Files.deleteIfExists(Paths.get(unsortedBam));
            }
            catch(IOException e)
            {
                BT_LOGGER.error("error deleting interim file: {}", e.toString());
            }
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AltContigRemapperConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AltContigRemapperConfig config = new AltContigRemapperConfig(configBuilder);
        new AltContigRemapper(config).run();
    }
}
