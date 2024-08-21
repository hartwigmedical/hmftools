package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.BreakendMatcher;
import com.hartwig.hmftools.esvee.utils.vcfcompare.match.VcfType;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SvVcfComparer
{
    private final String mSampleId;
    private final String mReferenceId;

    private final String mOldVcf;
    private final String mNewVcf;

    private final String mOldUnfilteredVcf;
    private final String mNewUnfilteredVcf;

    private final String mOutputDir;
    private final String mOutputId;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMap;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMap;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMapUnfiltered;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMapUnfiltered;

    private final boolean mShowNonPass;

    private final BreakendMatcher mBreakendMatcher;

    private RefGenomeVersion mRefGenomeVersion;

    private static final String OLD_VCF = "old_vcf";
    private static final String NEW_VCF = "new_vcf";

    private static final String OLD_UNFILTERED_VCF = "old_unfiltered_vcf";
    private static final String NEW_UNFILTERED_VCF = "new_unfiltered_vcf";

    private static final String SHOW_NON_PASS = "show_non_pass";

    public SvVcfComparer(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mReferenceId = configBuilder.getValue(REFERENCE, "");

        mOldVcf = configBuilder.getValue(OLD_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);

        mOldUnfilteredVcf = configBuilder.getValue(OLD_UNFILTERED_VCF);
        mNewUnfilteredVcf = configBuilder.getValue(NEW_UNFILTERED_VCF);

        mShowNonPass = configBuilder.hasFlag(SHOW_NON_PASS);

        mOutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mOldChrBreakendMap = new HashMap<>();
        mNewChrBreakendMap = new HashMap<>();

        mOldChrBreakendMapUnfiltered = new HashMap<>();
        mNewChrBreakendMapUnfiltered = new HashMap<>();

        mRefGenomeVersion = RefGenomeVersion.V37; // FIXME: Make this configurable

        mBreakendMatcher = new BreakendMatcher(mSampleId, mOutputDir, mOutputId, mRefGenomeVersion, mShowNonPass);
    }

    public void run()
    {
        if(mOldVcf == null || mNewVcf == null)
        {
            SV_LOGGER.error("Missing VCFs");
            return;
        }

        loadVariants(mOldVcf, mOldChrBreakendMap);
        loadVariants(mNewVcf, mNewChrBreakendMap);

        mBreakendMatcher.matchBreakends(mOldChrBreakendMap, mNewChrBreakendMap);

        if(mOldUnfilteredVcf != null)
        {
            loadVariants(mOldUnfilteredVcf, mOldChrBreakendMapUnfiltered);
            mBreakendMatcher.matchBreakends(mNewChrBreakendMap, mOldChrBreakendMapUnfiltered);
        }

        if(mNewUnfilteredVcf != null)
        {
            loadVariants(mNewUnfilteredVcf, mNewChrBreakendMapUnfiltered);
            mBreakendMatcher.matchBreakends(mOldChrBreakendMap, mNewChrBreakendMapUnfiltered);
        }

        mBreakendMatcher.gatherUnmatchedVariants(mOldChrBreakendMap, true);
        mBreakendMatcher.gatherUnmatchedVariants(mNewChrBreakendMap, false);
        mBreakendMatcher.writeBreakends();
        mBreakendMatcher.closeWriter();

        SV_LOGGER.info("Esvee compare VCFs complete");
    }

    private void loadVariants(final String vcfFile, final Map<String,List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.info("Loading vcfFile({})", vcfFile);

        VcfFileReader reader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = GenotypeIds.fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        SV_LOGGER.info("Genotype info: ref({}: {}) tumor({}: {})",
                genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);

        String currentChr = "";
        List<VariantBreakend> breakends = null;

        SvCaller svCaller = SvCaller.fromVcfPath(vcfFile);
        VcfType sourceVcfType = VcfType.fromVcfPath(vcfFile);

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

//            if(mRefGenomeVersion == null)
//                mRefGenomeVersion = chromosome.startsWith(CHR_PREFIX) ? V38 : V37;

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = new ArrayList<>();
                chrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(new VariantBreakend(variantContext, svCaller, sourceVcfType));
        }

        SV_LOGGER.info("Loaded {} SVs from {})",
                chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(), vcfFile);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addPath(OLD_VCF, true, "Path to the old VCF file");
        configBuilder.addPath(NEW_VCF, true, "Path to the new VCF file");

        configBuilder.addPath(OLD_UNFILTERED_VCF, false, "Path to the old unfiltered VCF file");
        configBuilder.addPath(NEW_UNFILTERED_VCF, false, "Path to the new unfiltered VCF file");

        configBuilder.addPath(SHOW_NON_PASS, false, "Show variants not PASSing in both old nor new VCF files");

        FileWriterUtils.addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        ConfigUtils.setLogLevel(configBuilder);

        SvVcfComparer svVcfCompare = new SvVcfComparer(configBuilder);
        svVcfCompare.run();
    }
}