package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import htsjdk.samtools.*;
import htsjdk.samtools.util.StringUtil;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.concurrent.atomic.AtomicInteger;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

public class Remapper
{
    private final RemapperConfig mConfig;

    public Remapper(final ConfigBuilder configBuilder)
    {
        mConfig = new RemapperConfig(configBuilder);
    }

    static boolean isRelevantDictionaryItem(SAMSequenceRecord samSequenceRecord) {
        return !samSequenceRecord.getSequenceName().toLowerCase().startsWith("hla");
    }

    static boolean hasAltReference(SAMRecord record) {
        return record.getReferenceName().toLowerCase().startsWith("hla") || record.getMateReferenceName().toLowerCase().startsWith("hla");
    }

    public void run()
    {
        BT_LOGGER.info("starting Remapper");

        long startTimeMs = System.currentTimeMillis();


//        BwaAligner aligner = new BwaAligner(mConfig.RefGenomeFile);

        try(SamReader samReader = SamReaderFactory. makeDefault().open(new File(mConfig.OrigBamFile))) {
            SAMFileHeader fileHeader = samReader.getFileHeader();
            BT_LOGGER.info("Header: " + fileHeader);
            SAMFileHeader unalteredRecordsHeader = fileHeader.clone();

            SAMSequenceDictionary dictionaryWithAltsRemoved = new SAMSequenceDictionary();
            fileHeader.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(Remapper::isRelevantDictionaryItem)
                    .forEach(dictionaryWithAltsRemoved::addSequence);
            unalteredRecordsHeader.setSequenceDictionary(dictionaryWithAltsRemoved);
//            unalteredRecordsHeader.getSequenceDictionary().
//            try (SAMFileWriter unalteredRecordsWriter = new SAMFileWriterFactory().makeBAMWriter(unalteredRecordsHeader, false, new File(mConfig.OutputFile))) {
            SAMFileWriter unalteredRecordsWriter = new SAMFileWriterFactory().makeBAMWriter(unalteredRecordsHeader, false, new File(mConfig.OutputFile));

//            long recordCount = samReader.iterator().stream().count();
//            BT_LOGGER.info("Total records: " + recordCount);
                AtomicInteger unalteredRecordCount = new AtomicInteger();
                samReader.forEach(record -> {
                    String bases = StringUtil.bytesToString(record.getReadBases());
//                if (bases.equals("GGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACCGCCTCTGCGGGGAGAAGCAAGGGGCCCTCCTGGCGGGGGCGCAGGACCGGGGGAGCCGCGCCGGGAGGAGGGTCGGGCAGGTCTCAGCCACTGCTCGCCCCC")) {
                    if (hasAltReference(record)) {
                    /*
                    BT_LOGGER.info("record: " + record);
                    List<BwaMemAlignment> alignments = aligner.alignSequence(record.getReadBases());
                    if (alignments.isEmpty()) {
//                        BT_LOGGER.info("no alignment found for " + record);
                    } else if (alignments.size() > 1) {
                        BT_LOGGER.info("multiple alignments found for " + record);

                    } else {
                        BT_LOGGER.info("alignment found for " + record);
                    }

                     */
                    } else {
                        // Write it out unaltered.
                        try {
                            unalteredRecordsWriter.addAlignment(record);
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                        unalteredRecordCount.getAndIncrement();
                    }
                });
                unalteredRecordsWriter.close();
                BT_LOGGER.info("Number of unaltered records: " + unalteredRecordCount.get());
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
        RemapperConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        Remapper regionSlicer = new Remapper(configBuilder);
        regionSlicer.run();
    }
}
