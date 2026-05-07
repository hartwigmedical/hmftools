package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SpliceLiftBack
{
    private final SpliceLiftBackConfig mConfig;

    public SpliceLiftBack(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceLiftBackConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        // TODO: can make this better, just getting it to work for now (PERFORMANCE IMPROVEMENTS)
        List<ContigEntry> contigEntries = ContigSidecar.read(mConfig.ContigSidecarFile);
        Map<String, ContigEntry> contigByName = new HashMap<>();
        for(ContigEntry entry : contigEntries)
            contigByName.put(entry.contigName(), entry);

        RD_LOGGER.info("loaded {} contig entries", contigEntries.size());

        processBam(mConfig.InputBam, mConfig.RefGenomeFile, contigByName, mConfig.formOutputBam());

        RD_LOGGER.info("SpliceLiftBack complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    static void processBam(
            final String inputBam, final String refGenomeFile, final Map<String, ContigEntry> contigByName, final String outputBam)
    {
        int total = 0;
        int passThrough = 0;
        int lifted = 0;
        int droppedUnliftable = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(refGenomeFile))
                .open(new File(inputBam)))
        {
            SAMFileHeader header = samReader.getFileHeader();

            try(SAMFileWriter samWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(outputBam)))
            {
                SAMRecordIterator iter = samReader.iterator();
                while(iter.hasNext())
                {
                    SAMRecord read = iter.next();
                    ++total;

                    if(read.getReadUnmappedFlag())
                    {
                        samWriter.addAlignment(read);
                        ++passThrough;
                        continue;
                    }

                    ContigEntry entry = contigByName.get(read.getReferenceName());
                    if(entry == null)
                    {
                        samWriter.addAlignment(read);
                        ++passThrough;
                        continue;
                    }

                    ContigTranslator.TranslationResult result = ContigTranslator.translate(
                            entry, read.getAlignmentStart(), read.getCigar());

                    if(result == null)
                    {
                        ++droppedUnliftable;
                        continue;
                    }

                    read.setReferenceName(result.Chromosome);
                    read.setAlignmentStart(result.GenomicStart);
                    read.setCigar(result.GenomicCigar);

                    samWriter.addAlignment(read);
                    ++lifted;
                }
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("BAM I/O failed: {}", e.toString());
            throw new RuntimeException(e);
        }

        RD_LOGGER.info("processed reads: total({}) passThrough({}) lifted({}) droppedUnliftable({})",
                total, passThrough, lifted, droppedUnliftable);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SpliceLiftBack liftBack = new SpliceLiftBack(configBuilder);
        liftBack.run();
    }
}
