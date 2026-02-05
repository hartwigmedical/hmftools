package com.hartwig.hmftools.pave.pon_gen;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.perf.StringCache;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.pon.PonChrCache;
import com.hartwig.hmftools.common.variant.pon.PonVariantData;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.ClinvarChrCache;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class RegionPonTask
{
    private final ChrBaseRegion mRegion;
    private final PonConfig mConfig;
    private final List<String> mSampleVcfs;

    private final VariantCache mVariantCache;
    private final ClinvarChrCache mClinvarData;
    private final PonAnnotation mExistingPon;
    private final HotspotRegionCache mSomaticHotspots;
    private final HotspotRegionCache mGermlineHotspots;
    private final TranscriptRegionCache mTranscriptChrCache;

    public RegionPonTask(
            final PonConfig config, final ChrBaseRegion region, final List<String> sampleVcfs, final PonAnnotation existingPon,
            final ClinvarAnnotation clinvarAnnotation, final HotspotCache hotspotCache, final EnsemblDataCache ensemblDataCache)
    {
        mConfig = config;
        mRegion = region;
        mSampleVcfs = sampleVcfs;
        mExistingPon = existingPon;

        mVariantCache = new VariantCache();

        // extract annotation information just for this region
        mClinvarData = new ClinvarChrCache(mRegion.Chromosome, new StringCache());

        ClinvarChrCache chrClinvarData = clinvarAnnotation.getChromosomeCache(mRegion.Chromosome);

        if(chrClinvarData != null)
        {
            chrClinvarData.entries().stream().filter(x -> mRegion.containsPosition(x.Position)).forEach(x -> mClinvarData.addEntry(x));
        }

        List<SimpleVariant> somaticHotspots = hotspotCache.getChromosomeSomaticHotspots(mRegion.Chromosome);

        if(somaticHotspots != null)
        {
            mSomaticHotspots = new HotspotRegionCache(
                    somaticHotspots.stream().filter(x -> mRegion.containsPosition(x.position())).collect(Collectors.toList()));
        }
        else
        {
            mSomaticHotspots = new HotspotRegionCache(Collections.emptyList());
        }

        List<SimpleVariant> germlineHotspots = hotspotCache.getChromosomeGermlineHotspots(mRegion.Chromosome);

        if(germlineHotspots != null)
        {
            mGermlineHotspots = new HotspotRegionCache(
                    germlineHotspots.stream().filter(x -> mRegion.containsPosition(x.position())).collect(Collectors.toList()));
        }
        else
        {
            mGermlineHotspots = new HotspotRegionCache(Collections.emptyList());
        }

        mTranscriptChrCache = TranscriptRegionCache.from(ensemblDataCache, mRegion);

        PV_LOGGER.debug("region({}) initialised: clinvarEntries({}) hotspots(somatic={} germline={}) exonicRegions({})",
                mRegion, mClinvarData.entryCount(), mSomaticHotspots.entryCount(), mGermlineHotspots.entryCount(),
                mTranscriptChrCache.entryCount());
    }

    public ChrBaseRegion region() { return mRegion; }
    public List<VariantPonData> variants() { return mVariantCache.variants(); }

    private static final int SAMPLE_LOG_COUNT = 100;
    private static final int SAMPLE_VARIANT_LOG_COUNT = 1_000_000;
    private static final int VARIANT_LOG_COUNT = 1_000_000;

    public void run()
    {
        for(int i = 0; i < mSampleVcfs.size(); ++i)
        {
            String sampleVcf = mSampleVcfs.get(i);
            loadSampleVariants(sampleVcf);

            if(i > 0 && (i % SAMPLE_LOG_COUNT) == 0)
            {
                PV_LOGGER.trace("region({}) processed {} samples", mRegion, i);
            }
        }

        int cachedVariantCount = mVariantCache.variantCount();

        addExistingPonEntries();

        // filter & finalise variants
        filterVariants();

        PV_LOGGER.debug("region({}) complete, variants(filter={} cached={}})", mRegion, mVariantCache.variantCount(), cachedVariantCount);
    }

    private void resetSearchIndices()
    {
        mVariantCache.resetSearch();

        if(mClinvarData != null)
            mClinvarData.resetSearch();

        mSomaticHotspots.resetSearch();
        mGermlineHotspots.resetSearch();
        mTranscriptChrCache.resetSearch();
    }

    private void loadSampleVariants(final String sampleVcf)
    {
        resetSearchIndices();

        VcfFileReader vcfReader = new VcfFileReader(sampleVcf);

        int varCount = 0;
        int varFilteredCount = 0;

        int lastPosition = 0;
        String lastChromosome = "";
        String lastRef = "";
        String lastAlt = "";

        for(VariantContext variantContext : vcfReader.regionIterator(mRegion))
        {
            ++varCount;

            double qual = variantContext.getPhredScaledQual();

            if(qual < mConfig.QualCutoff)
                continue;

            double mqf = variantContext.getAttributeAsDouble(MAP_QUAL_FACTOR, 0);

            if(mqf < mConfig.MqfCutoff)
                continue;

            int refSampleAd = 0;

            if(mConfig.RefSampleGenoptypeIndex >= 0)
            {
                Genotype refGenotype = variantContext.getGenotype(mConfig.RefSampleGenoptypeIndex);
                refSampleAd = refGenotype.getAD()[1];
            }

            int position = variantContext.getStart();
            String chromosome = variantContext.getContig();
            String ref = variantContext.getReference().getBaseString();
            String alt = variantContext.getAlternateAlleles().get(0).toString();

            // ignore duplicates (eg with different read-contexts)
            if(chromosome.equals(lastChromosome) && lastPosition == position && lastRef.equals(ref) && lastAlt.equals(alt))
                continue;

            VariantPonData variantPonData = getOrCreateVariant(variantContext, chromosome, position, ref, alt);
            variantPonData.addSampleData(refSampleAd);

            ++varFilteredCount;

            lastChromosome = chromosome;
            lastPosition = position;
            lastRef = ref;
            lastAlt = alt;

            if(PV_LOGGER.isTraceEnabled() && (varCount % SAMPLE_VARIANT_LOG_COUNT) == 0)
            {
                PV_LOGGER.trace("sampleVcf({}) loaded {} variants", sampleVcf, varCount);
            }
        }

        PV_LOGGER.trace("region({}) sampleVcf({}) loaded variants({} excluded={})",
                mRegion, sampleVcf, varCount, varCount- varFilteredCount);
    }

    private VariantPonData getOrCreateVariant(
            final VariantContext variantContext, final String chromosome, final int position, final String ref, final String alt)
    {
        VariantPonData variant = mVariantCache.getOrCreateVariant(chromosome, position, ref, alt);

        if(variant.sampleCount() > 0)
            return variant;

        if((mVariantCache.variantCount() % VARIANT_LOG_COUNT) == 0)
        {
            PV_LOGGER.debug("region({}) cached {} variants", mRegion, mVariantCache.variantCount());
        }

        // annotate variant
        if(mClinvarData != null)
        {
            Pathogenicity pathogenicity = mClinvarData.findPathogenicity(variant);
            variant.setClinvarPathogenicity(pathogenicity != null ? pathogenicity : Pathogenicity.UNKNOWN );
        }

        if(mSomaticHotspots.matchesHotspot(variant.Position, variant.Ref, variant.Alt))
        {
            variant.markSomaticHotspot();
        }

        if(mGermlineHotspots.matchesHotspot(variant.Position, variant.Ref, variant.Alt))
        {
            variant.markGermlineHotspot();
        }

        // check proximity to coding a region
        if(mTranscriptChrCache.inOrNearExonicRegion(variant.Position))
        {
            variant.markInCodingRegion();
        }

        if(variant.isIndel())
            variant.setRepeatCount(variantContext.getAttributeAsInt(SageVcfTags.REPEAT_COUNT, 0));

        return variant;
    }

    private void addExistingPonEntries()
    {
        if(!mExistingPon.enabled())
            return;

        mVariantCache.resetSearch();

        PonChrCache ponChrCache = mExistingPon.getChromosomeCache(mRegion.Chromosome);

        Map<Integer,List<PonVariantData>> positionMap = ponChrCache.positionMap();

        List<Integer> existingPositions = positionMap.keySet().stream()
                .filter(x -> mRegion.containsPosition(x)).sorted().toList();

        for(Integer existingPosition : existingPositions)
        {
            List<PonVariantData> existingVariants = positionMap.get(existingPosition);

            for(PonVariantData existingVariant : existingVariants)
            {
                VariantPonData variantPonData = mVariantCache.getOrCreateVariant(
                        mRegion.Chromosome, existingPosition, existingVariant.Ref, existingVariant.Alt);

                variantPonData.markInBasePonCache();

                if(variantPonData.sampleCount() > 0)
                {
                    variantPonData.markInMultiplePonCaches();

                    if(existingVariant.Samples > variantPonData.sampleCount())
                    {
                        variantPonData.setSampleCount(existingVariant.Samples);
                        variantPonData.setTotalReadCount(existingVariant.TotalSampleReads);
                    }

                    int maxSampleReadCount = max(variantPonData.maxSampleReadCount(), existingVariant.MaxSampleReads);
                    variantPonData.setMaxSampleReadCount(maxSampleReadCount);
                }
                else
                {
                    variantPonData.setSampleCount(existingVariant.Samples);
                    variantPonData.setMaxSampleReadCount(existingVariant.MaxSampleReads);
                    variantPonData.setTotalReadCount(existingVariant.TotalSampleReads);
                }
            }
        }
    }

    private void filterVariants()
    {
        List<VariantPonData> filteredVariants = mVariantCache.variants().stream().filter(x -> x.sampleCount() >= mConfig.MinSamples).collect(Collectors.toList());
        mVariantCache.replaceVariants(filteredVariants);
    }
}
