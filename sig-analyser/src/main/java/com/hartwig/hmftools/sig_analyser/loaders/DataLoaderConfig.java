package com.hartwig.hmftools.sig_analyser.loaders;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleListFile;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class DataLoaderConfig
{
    public final String OutputDir;
    public final String OutputFileId;
    public final List<String> SampleIds;
    public boolean ApplySampleQC;

    public final VariantFilters Filters;

    private static final String SAMPLE_IDS = "sample_ids";
    private static final String APPLY_SAMPLE_QC = "apply_sample_qc";
    private static final String SUBCLONAL_MIN = "subclonal_min";
    private static final String SUBCLONAL_MAX = "subclonal_max";
    private static final String PLOIDY_MAX = "ploidy_max";
    private static final String PLOIDY_MIN = "ploidy_min";

    public DataLoaderConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);

        ApplySampleQC = cmd.hasOption(APPLY_SAMPLE_QC);

        Filters = new VariantFilters(cmd);

        SampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_IDS))
        {
            final String sampleIdsStr = cmd.getOptionValue(SAMPLE_IDS);
            if(sampleIdsStr.contains(".csv"))
            {
                SampleIds.addAll(loadSampleListFile(sampleIdsStr));
            }
            else
            {
                SampleIds.addAll(Arrays.stream(sampleIdsStr.split(";")).collect(Collectors.toList()));
            }
        }
    }

    public void loadSampleIds(final DatabaseAccess dbAccess)
    {
        if(!SampleIds.isEmpty())
            return;

        if(ApplySampleQC)
            SampleIds.addAll(dbAccess.readPurpleSampleListPassingQC(MIN_SAMPLE_PURITY));
        else
            SampleIds.addAll(dbAccess.readSomaticVariantSampleList());
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_IDS, true, "Optional: restrict to sample list");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(APPLY_SAMPLE_QC, false, "Check sample QC status etc");
        options.addOption(SUBCLONAL_MAX, true, "Optional: subclonal max threshold");
        options.addOption(SUBCLONAL_MIN, true, "Optional: subclonal min threshold");
        options.addOption(PLOIDY_MAX, true, "Optional: ploidy max threshold");
        options.addOption(PLOIDY_MIN, true, "Optional: ploidy min threshold");
    }

}
