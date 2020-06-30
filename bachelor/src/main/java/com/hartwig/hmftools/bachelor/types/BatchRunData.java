package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BatchRunData
{
    private final List<String> mVcfFiles;
    private final List<String> mBamFiles;

    // config
    private static final String VCF_FILES = "batch_vcf_files";
    private static final String BAM_FILES = "batch_bam_files";

    public BatchRunData(final CommandLine cmd)
    {
        mVcfFiles = Lists.newArrayList();
        mBamFiles = Lists.newArrayList();

        /*
            try (final Stream<Path> stream = Files.walk(sampleDataPath, 1, FileVisitOption.FOLLOW_LINKS).parallel())
            {
                mSampleDataDirectories.addAll(stream.filter(p -> p.toFile().isDirectory())
                        .filter(p -> !p.equals(sampleDataPath))
                        .map(BatchRunData::new)
                        .collect(Collectors.toList()));
            }
            catch (Exception e)
            {
                BACH_LOGGER.error("failed find directories for batch run: {}", e.toString());
            }

            BACH_LOGGER.info("Found {} batch directories", mSampleDataDirectories.size());

         */

        loadFileList(cmd.getOptionValue(VCF_FILES), mVcfFiles);
        loadFileList(cmd.getOptionValue(BAM_FILES), mBamFiles);
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(VCF_FILES, true, "File with sample VCF file paths");
        options.addOption(BAM_FILES, true, "File with sample BAM file paths");
    }

    private void loadFileList(final String filename, final List<String> files)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String filePath = fileReader.readLine();

            while (filePath != null)
            {
                // check file exists before accept it
                if(!Files.exists(Paths.get(filePath)))
                {
                    BACH_LOGGER.error("skipping file({}) not found", filePath);
                    return;
                }

                files.add(filePath);

                filePath = fileReader.readLine();
            }

            BACH_LOGGER.info("loaded {} files from {}", files.size(), filename);

        }
        catch (IOException e)
        {
            BACH_LOGGER.error("failed to read sample list input CSV file({}): {}", filename, e.toString());
        }
    }
}
