package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.common.FileCommon;

import htsjdk.samtools.*;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
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
        try(SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.OrigBamFile)))
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
            String interimOutputFileName = mConfig.OutputFile + ".unsorted";
            File interimOutputFile = new File(interimOutputFileName);
            SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(newHeader, false, interimOutputFile);

            final BwaHlaRecordAligner aligner = new BwaHlaRecordAligner(mConfig.aligner(), newHeader, mConfig.RefGenVersion);
            HlaTransformer transformer = new HlaTransformer(aligner);
            samReader.forEach(record ->
                    transformer.process(record).forEach(bamWriter::addAlignment));

            BT_LOGGER.info("Finished processing input file. Number of HLA records processed: " + transformer.numberOfHlaRecordsProcessed());
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

            // Write the records to file.
            bamWriter.close();
            BT_LOGGER.info("BAM Writer closed.");

            // If the samtools path has been provided, sort the output. Else simply rename the unsorted file.
            if(mConfig.BamToolPath != null)
            {
                BT_LOGGER.info("Output file is to be sorted...");
                FileCommon.writeSortedBam(interimOutputFileName, mConfig.OutputFile, mConfig.BamToolPath, 1);
                BT_LOGGER.info("Sorting complete.");
            }
            else
            {
                File outputFile = new File(mConfig.OutputFile);
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
        BT_LOGGER.info("Remapping complete, mins({})", runTimeMinsStr(startTimeMs));
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
