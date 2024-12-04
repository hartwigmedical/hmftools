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
        BT_LOGGER.info("starting alt contig remapper");
        long startTimeMs = System.currentTimeMillis();

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

            HlaTransformer transformer = new HlaTransformer(new BwaHlaRecordAligner(mConfig.aligner()));
            samReader.forEach(record ->
                    transformer.process(record).forEach(bamWriter::addAlignment));

            // Deal with any unmatched reads.
            // Don't map these - log an error and write them out as they are
            List<SAMRecord> unmatched = transformer.unmatchedRecords();
            if (!unmatched.isEmpty()) {
                BT_LOGGER.warn("Some HLA contig records were unmatched. " + unmatched);
                unmatched.forEach(bamWriter::addAlignment);
            }

            // Write the records to file.
            bamWriter.close();

            // If the samtools path has been provided, sort the output. Else simply rename the unsorted file.
            if (mConfig.BamToolPath != null) {
                FileCommon.writeSortedBam(interimOutputFileName, mConfig.OutputFile, mConfig.BamToolPath, 1);
            } else {
                File outputFile = new File(mConfig.OutputFile);
                boolean renamed = interimOutputFile.renameTo(outputFile);
                if (!renamed) {
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
        AltContigRemapper regionSlicer = new AltContigRemapper(config);
        regionSlicer.run();
    }
}
