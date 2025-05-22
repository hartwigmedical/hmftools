package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.esvee.caller.VcfWriter.GERMLINE_VCF_ID;
import static com.hartwig.hmftools.esvee.caller.VcfWriter.SOMATIC_VCF_ID;
import static com.hartwig.hmftools.esvee.caller.VcfWriter.UNFILTERED_VCF_ID;

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class CompareConfig
{
    public final String OldVcf;
    public final String NewVcf;

    public final String OldUnfilteredVcf;
    public final String NewUnfilteredVcf;

    public final String OutputDir;
    public final String OutputId;

    public final boolean IncludeNonPass;
    public final boolean WriteMatches;

    public static final String OLD_VCF = "old_vcf";
    public static final String NEW_VCF = "new_vcf";

    public static final String OLD_UNFILTERED_VCF = "old_unfiltered_vcf";
    public static final String NEW_UNFILTERED_VCF = "new_unfiltered_vcf";

    private static final String INCLUDE_NON_PASS = "include_non_pass";
    private static final String WRITE_MATCHES = "write_matches";

    private static final String RUN_LINE_ROUTINE = "compare_line";

    public CompareConfig(final ConfigBuilder configBuilder)
    {
        OldVcf = configBuilder.getValue(OLD_VCF);
        NewVcf = configBuilder.getValue(NEW_VCF);

        OldUnfilteredVcf = getUnfilteredVcfFilename(configBuilder.getValue(OLD_UNFILTERED_VCF), OldVcf);
        NewUnfilteredVcf = getUnfilteredVcfFilename(configBuilder.getValue(NEW_UNFILTERED_VCF), NewVcf);

        OutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        IncludeNonPass = configBuilder.hasFlag(INCLUDE_NON_PASS);
        WriteMatches = configBuilder.hasFlag(WRITE_MATCHES);
    }

    private static String getUnfilteredVcfFilename(final String unfilteredVcfConfig, final String mainVcfConfig)
    {
        if(unfilteredVcfConfig != null)
            return unfilteredVcfConfig;

        if(mainVcfConfig.contains(SOMATIC_VCF_ID) || mainVcfConfig.contains(GERMLINE_VCF_ID))
        {
            String fileId = mainVcfConfig.contains(SOMATIC_VCF_ID) ? SOMATIC_VCF_ID : GERMLINE_VCF_ID;
            String unfilteredVcf = mainVcfConfig.replaceAll(fileId, UNFILTERED_VCF_ID);

            if(Files.exists(Paths.get(unfilteredVcf)))
                return unfilteredVcf;
        }

        return null;
    }

    /*
    public String formFilename(final String fileType)
    {
        String appStage = ESVEE_FILE_ID + FILE_NAME_DELIM + "compare" + FILE_NAME_DELIM + fileType;
        return formOutputFile(OutputDir, SampleId, appStage, "tsv", OutputId);
    }
    */

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(OLD_VCF, true, "Path to the old VCF file");
        configBuilder.addPath(NEW_VCF, true, "Path to the new VCF file");

        configBuilder.addPath(OLD_UNFILTERED_VCF, false, "Path to the old unfiltered VCF file");
        configBuilder.addPath(NEW_UNFILTERED_VCF, false, "Path to the new unfiltered VCF file");

        configBuilder.addFlag(INCLUDE_NON_PASS, "Show variants not PASSing in both old nor new VCF files");
        configBuilder.addFlag(WRITE_MATCHES, "Write exact matches to output");

        FileWriterUtils.addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
