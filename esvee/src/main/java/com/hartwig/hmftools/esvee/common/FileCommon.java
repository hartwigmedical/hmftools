package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BAM_FILE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_DISC_STATS_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public final class FileCommon
{
    public static final String APP_NAME = "Esvee";

    public static final String INPUT_VCF = "input_vcf";
    public static final String INPUT_VCF_DESC = "Input VCF";

    // file IDs and names
    public static final String ESVEE_FILE_ID = "esvee";
    public static final String PREP_FILE_ID = "esvee.prep";

    public static final String PREP_DIR = "esvee_prep_dir";
    public static final String PREP_DIR_DESC = "Esvee prep input directory";

    public static final String RAW_VCF_SUFFIX = "raw" + VCF_ZIP_EXTENSION;
    public static final String DEPTH_VCF_SUFFIX = "ref_depth" + VCF_ZIP_EXTENSION;

    // common config
    public static final String JUNCTION_FILE = "junction_file";
    public static final String JUNCTION_FILE_DESC = "Esvee Prep junction file, default is to match by sample name";

    public static final String KNOWN_HOTSPOT_FILE = "known_hotspot_file";

    public static final String FILE_NAME_DELIM = ".";

    public static final String REF_GENOME_IMAGE_EXTENSION = ".img";

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        SequencingType.registerConfig(configBuilder);
    }

    public static void setSequencingType(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SEQUENCING_TYPE_CFG))
        {
            SvConstants.Sequencing = SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));
        }
    }

    public static List<String> parseSampleBamLists(final ConfigBuilder configBuilder, final String configItem)
    {
        if(!configBuilder.hasValue(configItem))
            return Lists.newArrayList();

        return Arrays.stream(configBuilder.getValue(configItem).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
    }

    public static List<String> parseSampleIds(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
            return parseSampleBamLists(configBuilder, SAMPLE);

        List<String> sampleIds = Lists.newArrayList();
        sampleIds.addAll(parseSampleBamLists(configBuilder, TUMOR));
        sampleIds.addAll(parseSampleBamLists(configBuilder, REFERENCE));
        return sampleIds;
    }

    public static List<String> parseBamFiles(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(BAM_FILE))
            return parseSampleBamLists(configBuilder, BAM_FILE);

        List<String> sampleIds = Lists.newArrayList();
        sampleIds.addAll(parseSampleBamLists(configBuilder, TUMOR_BAM));
        sampleIds.addAll(parseSampleBamLists(configBuilder, REFERENCE_BAM));
        return sampleIds;
    }

    public static String formPrepInputFilename(
            final String outputDir, final String sampleId, final String fileId, @Nullable final String outputId)
    {
        // load from the output-id specific file if exists, otherwise revert to one without it
        String outputIdFilename = formOutputFile(outputDir, sampleId, PREP_FILE_ID, fileId, outputId);

        if(Files.exists(Paths.get(outputIdFilename)))
            return outputIdFilename;

        return formOutputFile(outputDir, sampleId, PREP_FILE_ID, fileId, null);
    }

    public static List<String> formPrepBamFilenames(final String outputDir, final List<String> sampleIds)
    {
        List<String> bamFiles = Lists.newArrayListWithExpectedSize(sampleIds.size());

        for(String sampleId : sampleIds)
        {
            String bamFile = format("%s%s.%s.bam", outputDir, sampleId, PREP_FILE_ID);

            if(!Files.exists(Paths.get(bamFile)))
                return Collections.emptyList();

            bamFiles.add(bamFile);
        }

        return bamFiles;
    }

    public static String formEsveeInputFilename(
            final String outputDir, final String sampleId, final String fileId, @Nullable final String outputId)
    {
        return formOutputFile(outputDir, sampleId, ESVEE_FILE_ID, fileId, outputId);
    }

    public static String formFragmentLengthDistFilename(final String outputDir, final String sampleId, final String outputId)
    {
        return formPrepInputFilename(outputDir, sampleId, PREP_FRAG_LENGTH_FILE_ID, outputId);
    }

    public static String formDiscordantStatsFilename(final String outputDir, final String sampleId, final String outputId)
    {
        return formPrepInputFilename(outputDir, sampleId, PREP_DISC_STATS_FILE_ID, outputId);
    }

    public static String formOutputFile(
            final String outputDir, final String sampleId, final String appStage, final String fileType, @Nullable final String outputId)
    {
        String filename = outputDir;

        filename += sampleId + FILE_NAME_DELIM + appStage;

        if(outputId != null)
            filename += FILE_NAME_DELIM + outputId;

        filename += FILE_NAME_DELIM + fileType;

        return filename;
    }


    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }

    public static void deleteInterimFile(final String filename)
    {
        try
        {
            Files.deleteIfExists(Paths.get(filename));
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error deleting interim file: {}", e.toString());
        }
    }

    public static void writeSortedBam(final String unsortedBam, final String sortedBam, final String bamToolPath, final int threads)
    {
        if(bamToolPath == null)
            return;

        SV_LOGGER.debug("writing sorted BAM: {}", sortedBam);

        BamToolName toolName = fromPath(bamToolPath);

        boolean success = BamOperations.sortBam(toolName, bamToolPath, unsortedBam, sortedBam, threads);

        if(success && toolName == BamToolName.SAMTOOLS)
        {
            success = BamOperations.indexBam(toolName, bamToolPath, sortedBam, threads);
        }

        if(success)
            deleteInterimFile(unsortedBam);
    }
}
