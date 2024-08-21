package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SvVcfComparer
{
    private final CompareConfig mConfig;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMap;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMap;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMapUnfiltered;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMapUnfiltered;

    private final BreakendMatcher mBreakendMatcher;

    public SvVcfComparer(final CompareConfig config)
    {
        mConfig = config;

        mOldChrBreakendMap = new HashMap<>();
        mNewChrBreakendMap = new HashMap<>();

        mOldChrBreakendMapUnfiltered = new HashMap<>();
        mNewChrBreakendMapUnfiltered = new HashMap<>();

        mBreakendMatcher = new BreakendMatcher(config);
    }

    public void run()
    {
        if(mConfig.OldVcf == null || mConfig.NewVcf == null)
        {
            SV_LOGGER.error("Missing VCFs");
            return;
        }

        loadVariants(mConfig.OldVcf, mOldChrBreakendMap);
        loadVariants(mConfig.NewVcf, mNewChrBreakendMap);

        mBreakendMatcher.matchBreakends(mOldChrBreakendMap, mNewChrBreakendMap);

        if(mConfig.OldUnfilteredVcf != null)
        {
            loadVariants(mConfig.OldUnfilteredVcf, mOldChrBreakendMapUnfiltered);
            mBreakendMatcher.matchBreakends(mNewChrBreakendMap, mOldChrBreakendMapUnfiltered);
        }

        if(mConfig.NewUnfilteredVcf != null)
        {
            loadVariants(mConfig.NewUnfilteredVcf, mNewChrBreakendMapUnfiltered);
            mBreakendMatcher.matchBreakends(mOldChrBreakendMap, mNewChrBreakendMapUnfiltered);
        }

        mBreakendMatcher.gatherUnmatchedVariants(mOldChrBreakendMap, mNewChrBreakendMap);

        mBreakendMatcher.writeBreakends();
        mBreakendMatcher.closeWriter();

        SV_LOGGER.info("Esvee compare VCFs complete");
    }

    private void loadVariants(final String vcfFile, final Map<String,List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.info("Loading vcfFile({})", vcfFile);

        VcfFileReader reader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = GenotypeIds.fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.SampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        String currentChr = "";
        List<VariantBreakend> breakends = null;

        SvCaller svCaller = SvCaller.fromVcfPath(vcfFile);
        VcfType sourceVcfType = VcfType.fromVcfPath(vcfFile);

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = new ArrayList<>();
                chrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(new VariantBreakend(variantContext, svCaller, sourceVcfType));
        }

        SV_LOGGER.info("  Loaded {} SVs from {})",
                chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(),
                vcfFile
        );
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