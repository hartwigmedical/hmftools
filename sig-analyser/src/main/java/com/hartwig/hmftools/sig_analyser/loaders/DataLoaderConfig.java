package com.hartwig.hmftools.sig_analyser.loaders;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleListFile;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DataLoaderConfig
{
    public final String OutputDir;
    public final String OutputFileId;
    public final List<String> SampleIds;
    public boolean ApplySampleQC;
    public final Double PloidyMin;
    public final Double PloidyMax;
    public final Double SubclonalLikelihoodMin;
    public final Double SubclonalLikelihoodMax;

    private static final String SAMPLE_IDS = "sample_ids";
    private static final String APPLY_SAMPLE_QC = "apply_sample_qc";
    private static final String SUBCLONAL_MIN = "subclonal_min";
    private static final String SUBCLONAL_MAX = "subclonal_max";
    private static final String PLOIDY_MAX = "ploidy_max";
    private static final String PLOIDY_MIN = "ploidy_min";

    private static final Logger LOGGER = LogManager.getLogger(DataLoaderConfig.class);

    public DataLoaderConfig(final CommandLine cmd)
    {
        OutputDir = cmd.getOptionValue(OUTPUT_DIR);
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);

        ApplySampleQC = cmd.hasOption(APPLY_SAMPLE_QC);

        SubclonalLikelihoodMin = initialiseDoubleValue(cmd, SUBCLONAL_MIN);
        SubclonalLikelihoodMax = initialiseDoubleValue(cmd, SUBCLONAL_MAX);

        PloidyMin = initialiseDoubleValue(cmd, PLOIDY_MIN);
        PloidyMax = initialiseDoubleValue(cmd, PLOIDY_MAX);

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

    private static Double initialiseDoubleValue(final CommandLine cmd, final String config)
    {
        if(cmd.hasOption(config))
            return Double.parseDouble(cmd.getOptionValue(config, "0"));
        else
            return null;
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

    public boolean passesFilters(final SomaticVariant variant)
    {
        if(SubclonalLikelihoodMin != null && variant.subclonalLikelihood() < SubclonalLikelihoodMin)
            return false;

        if(SubclonalLikelihoodMax != null && variant.subclonalLikelihood() > SubclonalLikelihoodMax)
            return false;

        if(PloidyMin != null && variant.variantCopyNumber() < PloidyMin)
            return false;

        if(PloidyMax != null && variant.variantCopyNumber() > PloidyMax)
            return false;

        return true;
    }
}
