package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class PurityConfig
{
    public final List<String> PatientIds;
    public final List<String> SomaticVcfs;
    public final String PurpleDir;
    public final String CobaltDir;
    public final String OutputDir;
    public final String OutputId;

    private static final String PATIENT_IDS = "patient_ids";
    private static final String SOMATIC_VCFS = "somatic_vcfs";
    private static final String VCF_DIR = "vcf_dir";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String COBALT_DIR = "cobalt_dir";

    public PurityConfig(final CommandLine cmd)
    {
        PatientIds = Lists.newArrayList();
        SomaticVcfs = Lists.newArrayList();

        if(cmd.hasOption(PATIENT_IDS))
        {
            Arrays.stream(cmd.getOptionValue(PATIENT_IDS).split(",", -1)).forEach(x -> PatientIds.add(x));
        }

        String vcfDir = cmd.hasOption(VCF_DIR) ? checkAddDirSeparator(cmd.getOptionValue(VCF_DIR)) : null;

        if(cmd.hasOption(SOMATIC_VCFS))
        {
            String[] vcfs = cmd.getOptionValue(SOMATIC_VCFS).split(",", -1);

            for(String vcf : vcfs)
            {
                vcf = vcf.trim();

                if(vcfDir != null)
                    SomaticVcfs.add(vcfDir + vcf);
                else
                    SomaticVcfs.add(vcf);
            }
        }

        PurpleDir = checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR));
        CobaltDir = checkAddDirSeparator(cmd.getOptionValue(COBALT_DIR));
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(PATIENT_IDS, true, "Patient IDs, separated by ','");
        options.addOption(SOMATIC_VCFS, true, "Somatic VCF files, separated by ','");
        options.addOption(VCF_DIR, true, "Directory for somatic VCF files");
        options.addOption(PURPLE_DIR, true, "Sample Purple directory");
        options.addOption(COBALT_DIR, true, "Sample Cobalt directory");
        addOutputOptions(options);
        addLoggingOptions(options);
    }
}
