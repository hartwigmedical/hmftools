package com.hartwig.hmftools.esvee.vcfcompare.match;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.caller.VcfWriter.SOMATIC_VCF_ID;
import static com.hartwig.hmftools.esvee.caller.VcfWriter.UNFILTERED_VCF_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.NEW_UNFILTERED_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.NEW_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.OLD_UNFILTERED_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.OLD_VCF;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.vcfcompare.CompareConfig;

import org.jetbrains.annotations.NotNull;

public class BreakendMatchTask implements Runnable
{
    private final CompareConfig mConfig;
    private final BreakendMatcher mBreakendMatcher;

    public BreakendMatchTask(final CompareConfig config)
    {
        mConfig = config;

        mBreakendMatcher = new BreakendMatcher(config);
    }

    @Override
    public void run()
    {
        SV_LOGGER.info("running breakend comparison");

        Map<String,List<VariantBreakend>> oldChrBreakendMap = loadVariants(mConfig.OldVcf, OLD_VCF);
        Map<String,List<VariantBreakend>> newChrBreakendMap = loadVariants(mConfig.NewVcf, NEW_VCF);

        mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMap);

        String oldUnfilteredVcf = getUnfilteredVcfFilename(mConfig.OldUnfilteredVcf, mConfig.OldVcf);
        String newUnfilteredVcf = getUnfilteredVcfFilename(mConfig.NewUnfilteredVcf, mConfig.NewVcf);

        if(oldUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> oldChrBreakendMapUnfiltered = loadVariants(oldUnfilteredVcf, OLD_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(newChrBreakendMap, oldChrBreakendMapUnfiltered);
        }

        if(newUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> newChrBreakendMapUnfiltered = loadVariants(newUnfilteredVcf, NEW_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMapUnfiltered);
        }

        mBreakendMatcher.gatherUnmatchedVariants(oldChrBreakendMap, newChrBreakendMap);

        BreakendMatchWriter writer = new BreakendMatchWriter(mBreakendMatcher.getBreakendMatches(), mConfig);
        writer.writeBreakends();

        SV_LOGGER.info("completed breakend comparison");
    }

    private String getUnfilteredVcfFilename(final String unfilteredVcfConfig, final String mainVcfConfig)
    {
        if(unfilteredVcfConfig != null)
            return unfilteredVcfConfig;

        String unfilteredVcf = mainVcfConfig.replaceAll(SOMATIC_VCF_ID, UNFILTERED_VCF_ID);

        return Files.exists(Paths.get(unfilteredVcf)) ? unfilteredVcf : null;
    }

    private static Map<String, List<VariantBreakend>> loadVariants(final String vcfFile, final String label)
    {
        return new SvVcfFile(vcfFile, label.toUpperCase()).loadVariants().getVariantsAsMap();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CompareConfig compareConfig = new CompareConfig(configBuilder);

        BreakendMatchTask svVcfCompare = new BreakendMatchTask(compareConfig);
        svVcfCompare.run();
    }
}