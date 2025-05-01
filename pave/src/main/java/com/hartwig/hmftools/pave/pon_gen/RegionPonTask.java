package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.ClinvarChrCache;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class RegionPonTask
{
    private final ChrBaseRegion mRegion;
    private final PonConfig mConfig;
    private final List<String> mSampleVcfs;

    private final ClinvarChrCache mClinvarData;
    private final HotspotRegionCache mSomaticHotspots;
    private final HotspotRegionCache mGermlineHotspots;
    private final TranscriptRegionCache mTranscriptChrCache;

    private final List<VariantPonData> mVariants;
    private int mLastVariantIndex; // to speed up searching

    public RegionPonTask(
            final PonConfig config, final ChrBaseRegion region, final List<String> sampleVcfs,
            final ClinvarAnnotation clinvarAnnotation, final HotspotCache hotspotCache, final EnsemblDataCache ensemblDataCache)
    {
        mConfig = config;
        mRegion = region;
        mSampleVcfs = sampleVcfs;

        // extract annotation information just for this region
        mClinvarData = new ClinvarChrCache(mRegion.Chromosome, new StringCache());

        ClinvarChrCache chrClinvarData = clinvarAnnotation.getChromosomeCache(mRegion.Chromosome);

        if(chrClinvarData != null)
        {
            chrClinvarData.entries().stream().filter(x -> mRegion.containsPosition(x.Position)).forEach(x -> mClinvarData.addEntry(x));
        }

        List<VariantHotspot> somaticHotspots = hotspotCache.getChromosomeSomaticHotspots(mRegion.Chromosome);

        if(somaticHotspots != null)
        {
            mSomaticHotspots = new HotspotRegionCache(
                    somaticHotspots.stream().filter(x -> mRegion.containsPosition(x.position())).collect(Collectors.toList()));
        }
        else
        {
            mSomaticHotspots = new HotspotRegionCache(Collections.emptyList());
        }

        List<VariantHotspot> germlineHotspots = hotspotCache.getChromosomeGermlineHotspots(mRegion.Chromosome);

        if(germlineHotspots != null)
        {
            mGermlineHotspots = new HotspotRegionCache(
                    germlineHotspots.stream().filter(x -> mRegion.containsPosition(x.position())).collect(Collectors.toList()));
        }
        else
        {
            mGermlineHotspots = new HotspotRegionCache(Collections.emptyList());
        }

        mTranscriptChrCache = new TranscriptRegionCache(ensemblDataCache, mRegion);

        PV_LOGGER.debug("region({}) initialised: clinvarEntries({}) hotspots(somatic={} germline={}) exonicRegions({})",
                mRegion, mClinvarData.entryCount(), mSomaticHotspots.entryCount(), mGermlineHotspots.entryCount(),
                mTranscriptChrCache.entryCount());

        mVariants = Lists.newArrayList();
        mLastVariantIndex = 0;
    }

    public ChrBaseRegion region() { return mRegion; }
    public List<VariantPonData> variants() { return mVariants; }

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

        PV_LOGGER.debug("region({}) VCFs loaded and cached variants({})", mRegion, mVariants.size());

        // filter & finalise variants
        filterVariants();

        PV_LOGGER.debug("region({}) complete, filtered variants({})", mRegion, mVariants.size());
    }

    private void loadSampleVariants(final String sampleVcf)
    {
        if(mClinvarData != null)
            mClinvarData.resetSearch();

        mSomaticHotspots.resetSearch();
        mGermlineHotspots.resetSearch();
        mTranscriptChrCache.resetSearch();

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

            if(!mConfig.ApplyWgsFilters)
            {
                String tier = variantContext.getAttributeAsString(TIER, "");
                if(tier.equals(HOTSPOT.toString()))
                    continue;
            }

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
        // start from the last inserted index since each VCF is ordered
        int index = mLastVariantIndex;
        while(index < mVariants.size())
        {
            VariantPonData variant = mVariants.get(index);

            if(position > variant.Position)
            {
                ++index;
                continue;
            }

            if(position < variant.Position)
                break;

            if(variant.Ref.equals(ref) && variant.Alt.equals(alt))
                return variant;

            ++index;
        }

        VariantPonData newVariant = new VariantPonData(chromosome, position, ref, alt);
        mVariants.add(index, newVariant);

        if((mVariants.size() % VARIANT_LOG_COUNT) == 0)
        {
            PV_LOGGER.debug("region({}) cached {} variants", mRegion, mVariants.size());
        }

        // annotate variant
        if(mClinvarData != null)
        {
            Pathogenicity pathogenicity = mClinvarData.findPathogenicity(newVariant);
            newVariant.setClinvarPathogenicity(pathogenicity != null ? pathogenicity : Pathogenicity.UNKNOWN );
        }

        if(mSomaticHotspots.matchesHotspot(newVariant.Position, newVariant.Ref, newVariant.Alt))
        {
            newVariant.markSomaticHotspot();
        }

        if(mGermlineHotspots.matchesHotspot(newVariant.Position, newVariant.Ref, newVariant.Alt))
        {
            newVariant.markGermlineHotspot();
        }

        // check proximity to coding a region
        if(mTranscriptChrCache.inOrNearExonicRegion(newVariant.Position))
        {
            newVariant.markInCodingRegion();
        }

        if(newVariant.isIndel())
            newVariant.setRepeatCount(variantContext.getAttributeAsInt(SageVcfTags.REPEAT_COUNT, 0));

        mLastVariantIndex = index;

        return newVariant;
    }

    private void filterVariants()
    {
        List<VariantPonData> filteredVariants = mVariants.stream().filter(x -> x.sampleCount() >= mConfig.MinSamples).collect(Collectors.toList());
        mVariants.clear();
        mVariants.addAll(filteredVariants);

        /*
        // only min sample count for now
        int index = 0;

        while(index < mVariants.size())
        {
            VariantPonData variant = mVariants.get(index);

            if(variant.sampleCount() < mConfig.MinSamples)
                mVariants.remove(index);
            else
                ++index;
        }
        */
    }
}
