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
import com.hartwig.hmftools.esvee.vcfcompare.line.LineLinkType;
import com.hartwig.hmftools.esvee.vcfcompare.line.LineLinkWriter;
import com.hartwig.hmftools.esvee.vcfcompare.line.LineLinker;
import com.hartwig.hmftools.esvee.vcfcompare.match.BreakendMatcher;

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
        SV_LOGGER.info("running LINE comparison");

        Map<String,List<VariantBreakend>> oldChrBreakendMap = loadAndLinkVariants(mConfig.OldVcf, OLD_VCF);
        Map<String,List<VariantBreakend>> newChrBreakendMap = loadAndLinkVariants(mConfig.NewVcf, NEW_VCF);

        mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMap, false);

        if(mConfig.OldUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> oldChrBreakendMapUnfiltered = loadAndLinkVariants(mConfig.OldUnfilteredVcf, OLD_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(newChrBreakendMap, oldChrBreakendMapUnfiltered, false);
        }

        if(mConfig.NewUnfilteredVcf != null)
        {
            Map<String,List<VariantBreakend>> newChrBreakendMapUnfiltered = loadAndLinkVariants(mConfig.NewUnfilteredVcf, NEW_UNFILTERED_VCF);
            mBreakendMatcher.matchBreakends(oldChrBreakendMap, newChrBreakendMapUnfiltered, false);
        }

        mBreakendMatcher.gatherUnmatchedVariants(oldChrBreakendMap, newChrBreakendMap);

        LineLinker.inferLinksBetweenBreakendSets(oldChrBreakendMap, newChrBreakendMap, LineLinkType.OLD_POLY_A_NEW_OTHER);
        LineLinker.inferLinksBetweenBreakendSets(newChrBreakendMap, oldChrBreakendMap, LineLinkType.NEW_POLY_A_OLD_OTHER);

        LineLinkWriter writer = new LineLinkWriter(mBreakendMatcher, mConfig);
        writer.writeBreakends();

        SV_LOGGER.info("completed LINE comparison");
    }

    private static Map<String, List<VariantBreakend>> loadAndLinkVariants(String vcfFile, String label)
    {
        Map<String, List<VariantBreakend>> chrBreakendMap =  new SvVcfFile(vcfFile, label.toUpperCase())
                .loadVariants()
                .getVariantsAsMap();

        Map<String, List<VariantBreakend>> chrBreakendMapDeduped = LineLinker.dedupBreakends(chrBreakendMap);

        LineLinker.linkBreakends(chrBreakendMapDeduped);

        return chrBreakendMapDeduped;
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
