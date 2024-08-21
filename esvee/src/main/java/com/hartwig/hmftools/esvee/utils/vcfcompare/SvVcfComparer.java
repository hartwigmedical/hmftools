package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.NotNull;

public class SvVcfComparer
{
    private final CompareConfig mConfig;
    private final BreakendMatcher mBreakendMatcher;

    public SvVcfComparer(final CompareConfig config)
    {
        mConfig = config;

        mBreakendMatcher = new BreakendMatcher(config);
    }

    public void run()
    {
        if(mConfig.OldVcf == null || mConfig.NewVcf == null)
        {
            SV_LOGGER.error("Missing VCFs");
            return;
        }

        Map<String,List<VariantBreakend>> oldChrBreakendMap = VariantBreakend.loadVariants(mConfig.OldVcf);
        Map<String,List<VariantBreakend>> newChrBreakendMap = VariantBreakend.loadVariants(mConfig.NewVcf);

        mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMap);

        if(mConfig.OldUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> oldChrBreakendMapUnfiltered = VariantBreakend.loadVariants(mConfig.OldUnfilteredVcf);
            mBreakendMatcher.matchBreakends(newChrBreakendMap, oldChrBreakendMapUnfiltered);
        }

        if(mConfig.NewUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> newChrBreakendMapUnfiltered = VariantBreakend.loadVariants(mConfig.NewUnfilteredVcf);
            mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMapUnfiltered);
        }

        mBreakendMatcher.gatherUnmatchedVariants(oldChrBreakendMap, newChrBreakendMap);

        mBreakendMatcher.writeBreakends();
        mBreakendMatcher.closeWriter();

        SV_LOGGER.info("Esvee compare VCFs complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CompareConfig compareConfig = new CompareConfig(configBuilder);

        SvVcfComparer svVcfCompare = new SvVcfComparer(compareConfig);
        svVcfCompare.run();
    }
}