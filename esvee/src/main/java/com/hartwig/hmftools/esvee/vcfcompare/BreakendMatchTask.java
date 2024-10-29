package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.NEW_UNFILTERED_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.NEW_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.OLD_UNFILTERED_VCF;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareConfig.OLD_VCF;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.vcfcompare.common.SvVcfFile;
import com.hartwig.hmftools.esvee.vcfcompare.common.VariantBreakend;
import com.hartwig.hmftools.esvee.vcfcompare.match.BreakendMatchWriter;
import com.hartwig.hmftools.esvee.vcfcompare.match.BreakendMatcher;

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
        SV_LOGGER.info("Running task: " + CompareTask.MATCH_BREAKENDS);

        Map<String,List<VariantBreakend>> oldChrBreakendMap = loadVariants(mConfig.OldVcf, OLD_VCF);
        Map<String,List<VariantBreakend>> newChrBreakendMap = loadVariants(mConfig.NewVcf, NEW_VCF);

        mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMap);

        if(mConfig.OldUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> oldChrBreakendMapUnfiltered = loadVariants(mConfig.OldUnfilteredVcf, OLD_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(newChrBreakendMap, oldChrBreakendMapUnfiltered);
        }

        if(mConfig.NewUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> newChrBreakendMapUnfiltered = loadVariants(mConfig.NewUnfilteredVcf, NEW_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMapUnfiltered);
        }

        mBreakendMatcher.gatherUnmatchedVariants(oldChrBreakendMap, newChrBreakendMap);

        BreakendMatchWriter writer = new BreakendMatchWriter(mBreakendMatcher.getBreakendMatches(), mConfig);
        writer.writeBreakends();

        SV_LOGGER.info("Completed task: " + CompareTask.MATCH_BREAKENDS);
    }

    private static Map<String, List<VariantBreakend>> loadVariants(String vcfFile, String label)
    {
        return new SvVcfFile(vcfFile, label.toUpperCase())
                .loadVariants()
                .getVariantsAsMap();
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