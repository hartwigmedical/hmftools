package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.utils.vcfcompare.line.LineLinkWriter;
import com.hartwig.hmftools.esvee.utils.vcfcompare.line.LineLinker;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.NotNull;

public class LineCompareTask implements Runnable
{
    private final CompareConfig mConfig;
    private final BreakendMatcher mBreakendMatcher;

    public LineCompareTask(final CompareConfig config)
    {
        mConfig = config;

        mBreakendMatcher = new BreakendMatcher(config);
    }

    @Override
    public void run()
    {
        SV_LOGGER.info("Running task: " + CompareTask.LINE_COMPARE);

        Map<String, List<VariantBreakend>> oldChrBreakendMap = loadAndLinkVariants(mConfig.OldVcf);
        Map<String, List<VariantBreakend>> newChrBreakendMap = loadAndLinkVariants(mConfig.NewVcf);

        mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMap, false);

        if(mConfig.OldUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> oldChrBreakendMapUnfiltered = loadAndLinkVariants(mConfig.OldUnfilteredVcf);
            mBreakendMatcher.matchBreakends(newChrBreakendMap, oldChrBreakendMapUnfiltered, false);
        }

        if(mConfig.NewUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> newChrBreakendMapUnfiltered = loadAndLinkVariants(mConfig.NewUnfilteredVcf);
            mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMapUnfiltered, false);
        }

        mBreakendMatcher.gatherUnmatchedVariants(oldChrBreakendMap, newChrBreakendMap);

        LineLinker.inferLinks(oldChrBreakendMap, newChrBreakendMap);
        LineLinker.inferLinks(newChrBreakendMap, oldChrBreakendMap);

        LineLinkWriter writer = new LineLinkWriter(mBreakendMatcher, mConfig);
        writer.writeBreakends();

        SV_LOGGER.info("Completed task: " + CompareTask.LINE_COMPARE);
    }

    public static Map<String, List<VariantBreakend>> loadAndLinkVariants(String vcfFile)
    {
        Map<String, List<VariantBreakend>> chrBreakendMap = VariantBreakend.loadVariants(vcfFile);

        LineLinker.linkBreakends(chrBreakendMap);

        return chrBreakendMap;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CompareConfig compareConfig = new CompareConfig(configBuilder);

        LineCompareTask svVcfCompare = new LineCompareTask(compareConfig);
        svVcfCompare.run();
    }
}
