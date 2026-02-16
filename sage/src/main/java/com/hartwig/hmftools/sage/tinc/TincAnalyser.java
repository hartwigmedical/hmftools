package com.hartwig.hmftools.sage.tinc;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.parseIntegerList;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_RELATIVE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilterConfig.getTieredSoftFilterConfig;
import static com.hartwig.hmftools.sage.filter.VariantFilters.aboveMaxGermlineRelativeQual;
import static com.hartwig.hmftools.sage.filter.VariantFilters.aboveMaxGermlineVaf;
import static com.hartwig.hmftools.sage.filter.VariantFilters.aboveMaxMnvIndelGermlineAltSupport;
import static com.hartwig.hmftools.sage.filter.VariantFilters.includeRefJitterInMsiIndel;
import static com.hartwig.hmftools.sage.filter.VariantFilters.isPanelOrLongIndelRepeatVariant;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.populateDefaultLevels;
import static com.hartwig.hmftools.sage.tinc.TincConstants.RECOVERY_FILTERS;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_RECOVERY_FACTOR;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_RECOVERY_MIN;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.SIMPLE_ALT_COUNT;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadSupportCounts;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.filter.SoftFilterConfig;

import com.hartwig.hmftools.sage.seqtech.UltimaUtils;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.logging.log4j.util.Strings;

public class TincAnalyser
{
    private final TincConfig mConfig;
    private final VariantCache mVariantCache;
    private double mCalculatedTincLevel;

    public TincAnalyser(final TincConfig config)
    {
        mConfig = config;
        mVariantCache = new VariantCache(mConfig);
        mCalculatedTincLevel = 0;
    }

    public void run(final FilterConfig filterConfig)
    {
        SG_LOGGER.info("running TINC analyser");

        mVariantCache.loadVariants();

        mCalculatedTincLevel = TincCalculator.calculate(mVariantCache.fittingVariants(), populateDefaultLevels());

        if(filterConfig != null && mCalculatedTincLevel > 0)
        {
            recoverVariants(filterConfig, mVariantCache.variants(), mCalculatedTincLevel);
        }
    }

    public static String generateTincVcfFilename(final String inputVcf)
    {
        return inputVcf.replace(VCF_ZIP_EXTENSION, ".tinc" + VCF_ZIP_EXTENSION);
    }

    public void writeVcf(final IndexedFastaSequenceFile refGenome, final String inputVcf, final String outputVcf)
    {
        TincVcfWriter.writeVcf(refGenome, inputVcf, outputVcf, mVariantCache.variants(), mCalculatedTincLevel);
    }

    protected static boolean filterOutVariant(final VariantData variant)
    {
        double reducedAltFrags = variant.calcReducedAltValue(variant.ReferenceAltFrags);

        if(reducedAltFrags <= DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT)
            return false;

        if(variant.Context.hasAttribute(LOCAL_PHASE_SET))
            return false;

        return true;
    }

    @VisibleForTesting
    public static void recoverVariants(final FilterConfig filterConfig, final List<VariantData> variants, final double tincLevel)
    {
        /*
        - If a TINC > 0% is found, calculate a recoveryTinc = 2.5 * TINC + 3%. For example if TINC = 10%, then recoveryTinc = 28%. If TINC = 2% then recoveryTinc = 8%

            The intuitive explanation for why this needs a static 3% term is because of the discreteness of AD at low TINCs.
            If the expected germlineAD for a subset of variants for a given TINC is 0.1 then we can expect most of these variants to have 0 AD,
            but the ones that need to be recovered will have at least 1, i.e. 10x the amount expected from the TINC.
            Thus backing out a fixed factor of TINC is not sufficient

        - We back out the level of germline AF implied by this recoveryTinc.
            eg if tumorAF was 30% we would reduce germlineAF by 0.28 * 30% = 8.4%

        - If the resultant germlineAF is below our 4% threshold, and the variant is recoverable, then it is recovered.
        We do the same thing with the 4% maxGermlineRelRawQual threshold
         */

        double recoveryTinc = TINC_RECOVERY_MIN + tincLevel * TINC_RECOVERY_FACTOR;

        int recoveredCount = 0;

        for(VariantData variant : variants)
        {
            if(variant.isPassing())
                continue;

            // remove germline filters from any variant with them, not just those with only germline filters (ie recoverable)
            if(variant.Context.getFilters().stream().noneMatch(x -> RECOVERY_FILTERS.contains(x)))
                continue;

            // skip checking if the variant is in a PON
            if(variant.ponFiltered())
                continue;

            variant.buildNewFilters();

            double afReduction = min(variant.tumorAf() * recoveryTinc, 1.0);
            variant.setReferenceAltFragReduction(afReduction);

            VariantTier tier = VariantTier.fromContext(variant.Context);

            SoftFilterConfig tierConfig = getTieredSoftFilterConfig(tier, filterConfig);

            // extract values used in more than one filter
            ReadSupportCounts refReadCounts = new ReadSupportCounts(parseIntegerList(variant.RefGenotype, READ_CONTEXT_COUNT));

            recheckMaxGermlineVaf(tierConfig, variant, refReadCounts);
            recheckMaxGermlineRelativeQual(variant,refReadCounts);
            recheckMaxGermlineAltSupport(variant, tier, refReadCounts);

            if(variant.recovered())
                ++recoveredCount;
        }

        if(recoveredCount > 0)
        {
            SG_LOGGER.debug("recovered {} variants", recoveredCount);
        }
    }

    private static void recheckMaxGermlineVaf(
            final SoftFilterConfig config, final VariantData variant, final ReadSupportCounts refReadCounts)
    {
        if(!variant.newFilters().contains(MAX_GERMLINE_VAF))
            return;

        double tumorVaf = variant.tumorAf();

        int simpleAltMatches = getGenotypeAttributeAsInt(variant.RefGenotype, SIMPLE_ALT_COUNT, 0);

        int altSupport = refReadCounts.altSupport();
        int adjustedRefAltCount = altSupport + simpleAltMatches;
        int rcRepeatCount = variant.Context.getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0);
        String rcRepeatBases = variant.Context.getAttributeAsString(READ_CONTEXT_REPEAT_SEQUENCE, Strings.EMPTY);

        if(variant.isLongIndel() || includeRefJitterInMsiIndel(rcRepeatCount, rcRepeatBases, variant.Ref, variant.Alt))
        {
            int[] jitterCounts = parseIntegerList(variant.RefGenotype, READ_CONTEXT_JITTER);
            adjustedRefAltCount += jitterCounts[0] + jitterCounts[1];
        }

        adjustedRefAltCount = variant.calcReducedAltCount(adjustedRefAltCount);

        double qual = variant.Context.getPhredScaledQual();
        boolean isValidIndelRepeat = isPanelOrLongIndelRepeatVariant(variant.tier(), qual, variant.isIndel(), rcRepeatCount);

        if(aboveMaxGermlineVaf(variant.tier(), isValidIndelRepeat, tumorVaf, adjustedRefAltCount, refReadCounts.Total, config.MaxGermlineVaf))
            return;

        variant.newFilters().remove(MAX_GERMLINE_VAF);
    }

    private static void recheckMaxGermlineRelativeQual(final VariantData variant, final ReadSupportCounts refReadCounts)
    {
        if(!variant.newFilters().contains(MAX_GERMLINE_RELATIVE_QUAL))
            return;

        int[] avgTumorBaseQuals = parseIntegerList(variant.TumorGenotype, AVG_RECALIBRATED_BASE_QUAL);
        int[] avgRefBaseQuals = parseIntegerList(variant.RefGenotype, AVG_RECALIBRATED_BASE_QUAL);

        ReadSupportCounts tumorReadQuals = new ReadSupportCounts(parseIntegerList(variant.TumorGenotype, READ_CONTEXT_QUALITY));
        ReadSupportCounts refReadQuals = new ReadSupportCounts(parseIntegerList(variant.RefGenotype, READ_CONTEXT_QUALITY));
        ReadSupportCounts tumorReadCounts = new ReadSupportCounts(parseIntegerList(variant.TumorGenotype, READ_CONTEXT_COUNT));

        double tumorVaf = variant.tumorAf();
        int tumorDepth = tumorReadCounts.Total;
        int tumorAvgBaseQual = avgTumorBaseQuals[1];
        int refDepth = refReadCounts.Total;
        int refAltSupport = refReadCounts.altSupport();
        int refAvgBaseQual = avgRefBaseQuals[1];

        double tumorQual = tumorReadQuals.Full + tumorReadQuals.PartialCore + tumorReadQuals.Realigned;
        double refQual = refReadQuals.Full + refReadQuals.PartialCore + refReadQuals.Realigned;

        refQual = variant.calcReducedAltValue(refQual);

        double qual = variant.Context.getPhredScaledQual();
        int rcRepeatCount = variant.Context.getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0);
        boolean isValidIndelRepeat = isPanelOrLongIndelRepeatVariant(variant.tier(), qual, variant.isIndel(), rcRepeatCount);

        if(aboveMaxGermlineRelativeQual(
                variant.tier(), tumorQual, tumorVaf, tumorDepth, tumorAvgBaseQual, refQual, refDepth, refAltSupport, refAvgBaseQual, isValidIndelRepeat))
            return;

        variant.newFilters().remove(MAX_GERMLINE_RELATIVE_QUAL);
    }

    private static void recheckMaxGermlineAltSupport(final VariantData variant, final VariantTier tier, final ReadSupportCounts refReadCounts)
    {
        if(!variant.newFilters().contains(MAX_GERMLINE_ALT_SUPPORT))
            return;

        boolean isLongInsert = variant.isInsert() && variant.isLongIndel();

        int newAltSupport = variant.calcReducedAltCount(refReadCounts.altSupport());

        if(aboveMaxMnvIndelGermlineAltSupport(tier, variant.isMnv(), isLongInsert, refReadCounts.Total, newAltSupport))
            return;

        variant.newFilters().remove(MAX_GERMLINE_ALT_SUPPORT);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        TincConfig.registerAppConfig(configBuilder);
        FilterConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        long startTimeMs = System.currentTimeMillis();

        TincConfig config = new TincConfig(configBuilder);
        TincAnalyser tincAnalyser = new TincAnalyser(config);

        FilterConfig filterConfig = new FilterConfig(configBuilder);
        tincAnalyser.run(filterConfig);

        if(config.RewriteVcf && configBuilder.hasValue(REF_GENOME))
        {
            String outputVcf = generateTincVcfFilename(config.InputVcf);

            SG_LOGGER.debug("writing TINC VCF: {}", outputVcf);
            tincAnalyser.writeVcf(loadRefGenome(configBuilder.getValue(REF_GENOME)).refGenomeFile(), config.InputVcf, outputVcf);
        }

        SG_LOGGER.info("TINC complete, mins({})", runTimeMinsStr(startTimeMs));
    }
}
