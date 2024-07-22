package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.PurpleConfig;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantCache
{
    private final PurpleConfig mConfig;

    private final List<SomaticVariant> mVariants;

    private VCFHeader mVcfHeader;
    private GenotypeIds mGenotypeIds;

    // counts for plot & chart down-sampling
    private int mIndelCount;
    private int mSnpCount;

    public SomaticVariantCache(final PurpleConfig config)
    {
        mConfig = config;

        mVariants = Lists.newArrayList();
        mIndelCount = 0;
        mSnpCount = 0;
        mVcfHeader = null;
        mGenotypeIds = null;
    }

    public boolean hasData() { return !mVariants.isEmpty(); }
    public List<SomaticVariant> variants() { return mVariants; }
    public GenotypeIds genotypeIds() { return mGenotypeIds; }

    public int snpCount() { return mSnpCount; }
    public int indelCount() { return mIndelCount; }

    public void loadSomatics(final String somaticVcf, final ListMultimap<Chromosome,VariantHotspot> somaticHotspots)
    {
        if(somaticVcf.isEmpty())
            return;

        final HotspotEnrichment hotspotEnrichment = new HotspotEnrichment(somaticHotspots, true);

        VcfFileReader vcfReader = new VcfFileReader(somaticVcf);
        mVcfHeader = vcfReader.vcfHeader();

        mGenotypeIds = GenotypeIds.fromVcfHeader(mVcfHeader, mConfig.ReferenceId, mConfig.TumorId);

        boolean tumorOnly = mConfig.tumorOnlyMode();

        for(VariantContext variantContext : vcfReader.iterator())
        {
            SomaticVariant variant = new SomaticVariant(variantContext, mConfig.TumorId, mConfig.ReferenceId);

            if(tumorOnly && HlaCommon.containsPosition(variant)) // ignore these completely
                continue;

            if(!mConfig.TierQualFilters.isEmpty())
            {
                Integer qualThreshold = mConfig.TierQualFilters.get(variant.decorator().tier());
                if(qualThreshold != null && variant.decorator().qual() < qualThreshold)
                    continue;
            }

            if(mConfig.FilterSomaticsOnGene)
            {
                if(variant.variantImpact() == null || variant.variantImpact().GeneName.isEmpty())
                    continue;
            }

            if(mConfig.excludeOnSpecificRegion(variant.chromosome(), variant.position()))
                continue;

            mVariants.add(variant);

            // hotspot status is used in fitting as well as during and for enrichment
            hotspotEnrichment.processVariant(variantContext);

            if(variant.isPass())
            {
                if(variant.type() == VariantType.INDEL)
                    mIndelCount++;
                else
                    mSnpCount++;
            }
        }

        PPL_LOGGER.info("loaded {} somatic variants from {}", mVariants.size(), somaticVcf);
    }

    public VCFHeader getVcfHeader() { return mVcfHeader; }

    public void purityEnrich(final SomaticPurityEnrichment purityEnrichment)
    {
        mVariants.forEach(x -> purityEnrichment.processVariant(x));
    }
}
