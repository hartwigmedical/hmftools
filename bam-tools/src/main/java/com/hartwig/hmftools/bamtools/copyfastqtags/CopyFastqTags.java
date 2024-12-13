package com.hartwig.hmftools.bamtools.copyfastqtags;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CopyFastqTags
{
    private final CopyFastqTagsConfig mConfig;
    private final Map<String, List<FastqTag>> mReadToTags;

    public CopyFastqTags(final ConfigBuilder configBuilder)
    {
        mConfig = new CopyFastqTagsConfig(configBuilder);
        mReadToTags = Maps.newHashMap();
    }

    public void run()
    {
        BT_LOGGER.info("starting CopyFastqTags, writing output to {}", mConfig.OutputBamFile);
        long startTimeMs = System.currentTimeMillis();

        processFastq();
        processBam();

        BT_LOGGER.info("BamCompare complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processFastq()
    {
        try(BufferedReader fastqReader = createBufferedReader(mConfig.FastqFile))
        {
            FastqReadTags fastqReadTags;
            while((fastqReadTags = nextFastqReadTags(fastqReader)) != null)
                mReadToTags.put(fastqReadTags.ReadName, fastqReadTags.Tags);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    private void processBam()
    {
        SamReader samReader = null;
        SAMFileWriter samWriter = null;
        try
        {
            samReader = SamReaderFactory.makeDefault().open(new File(mConfig.BamFile));
            SAMFileHeader header = samReader.getFileHeader();
            samWriter = new SAMFileWriterFactory().makeBAMWriter(header, true, new File(mConfig.OutputBamFile));

            SAMRecordIterator iter = samReader.iterator();
            while(iter.hasNext())
            {
                SAMRecord read = iter.next();
                if(read.getReadPairedFlag())
                    throw new RuntimeException(format("Read with readName(%s) is paired, only unpaired reads are supported", read.getReadName()));

                List<FastqTag> fastqTags = mReadToTags.get(read.getReadName());
                if(fastqTags == null)
                {
                    samWriter.addAlignment(read);
                    continue;
                }

                for(FastqTag tag : fastqTags)
                    read.setAttribute(tag.TagName, tag.Value);

                samWriter.addAlignment(read);
            }
        }
        finally
        {
            try
            {
                if(samReader != null)
                    samReader.close();

                if(samWriter != null)
                    samWriter.close();
            }
            catch(IOException e)
            {
                throw new RuntimeException(e);
            }
        }
    }

    @Nullable
    private static FastqReadTags nextFastqReadTags(final BufferedReader reader)
    {
        try
        {
            String readNameAndTags = reader.readLine();
            if(readNameAndTags == null)
                return null;

            String[] components = readNameAndTags.split("\\s+");
            String readName = components[0].substring(1);

            reader.readLine();
            reader.readLine();
            if(reader.readLine() == null)
                throw new RuntimeException(format("Partial fastq record found: readName(%s).", readName));

            List<FastqTag> tags = Lists.newArrayList();
            for(int i = 1; i < components.length; i++)
                tags.add(FastqTag.fromFormattedAttribute(components[i]));

            return new FastqReadTags(readName, tags);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CopyFastqTagsConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CopyFastqTags copyFastqTags = new CopyFastqTags(configBuilder);
        copyFastqTags.run();
    }
}
