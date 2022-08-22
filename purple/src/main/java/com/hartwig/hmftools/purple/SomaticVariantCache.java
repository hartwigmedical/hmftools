package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.somatic.HotspotEnrichment;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;
import com.hartwig.hmftools.purple.somatic.SomaticPurityEnrichment;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantCache
{
    private final PurpleConfig mConfig;

    private final CompoundFilter mFilter;

    private final List<SomaticVariant> mVariants;
    private final List<SomaticVariant> mFittingVariants; // pass and SNVs only

    private VCFHeader mVcfHeader;

    // counts for plot & chart down-sampling
    private int mIndelCount;
    private int mSnpCount;

    public SomaticVariantCache(final PurpleConfig config)
    {
        mConfig = config;

        mFilter = new CompoundFilter(true);
        mFilter.add(new SGTFilter());
        mFilter.add(new HumanChromosomeFilter());
        mFilter.add(new NTFilter());

        mVariants = Lists.newArrayList();
        mFittingVariants = Lists.newArrayList();
        mIndelCount = 0;
        mSnpCount = 0;
        mVcfHeader = null;
    }

    public boolean hasData() { return !mVariants.isEmpty(); }
    public List<SomaticVariant> variants() { return mVariants; }
    public List<SomaticVariant> fittingVariants() { return mFittingVariants; }

    public int snpCount() { return mSnpCount; }
    public int indelCount() { return mIndelCount; }

    public void loadSomatics(final String somaticVcf, final ListMultimap<Chromosome,VariantHotspot> somaticHotspots)
    {
        if(somaticVcf.isEmpty())
            return;

        final HotspotEnrichment hotspotEnrichment = new HotspotEnrichment(somaticHotspots, true);

        VCFFileReader vcfReader = new VCFFileReader(new File(somaticVcf), false);
        mVcfHeader = vcfReader.getHeader();

        boolean tumorOnly = mConfig.tumorOnlyMode();

        for(VariantContext variantContext : vcfReader)
        {
            SomaticVariant variant = new SomaticVariant(variantContext, mConfig.TumorId);

            if(tumorOnly && HlaCommon.containsPosition(variant)) // ignore these completely
                continue;

            if(!mConfig.TierQualFilters.isEmpty())
            {
                Integer qualThreshold = mConfig.TierQualFilters.get(variant.decorator().tier());
                if(qualThreshold != null && variant.decorator().qual() < qualThreshold)
                    continue;
            }

            mVariants.add(variant);

            // hotspot status is used in fitting as well as during and for enrichment
            hotspotEnrichment.processVariant(variantContext);

            if(!tumorOnly && isFittingCandidate(variant))
                mFittingVariants.add(variant);

            if(variant.isPass())
            {
                if(variant.type() == VariantType.INDEL)
                    mIndelCount++;
                else
                    mSnpCount++;
            }
        }

        PPL_LOGGER.info("loaded somatic variants({} fitting={}) from {}",
                mVariants.size(), mFittingVariants.size(), somaticVcf);
    }

    public VCFHeader getVcfHeader() { return mVcfHeader; }

    public void purityEnrich(final SomaticPurityEnrichment purityEnrichment)
    {
        mVariants.forEach(x -> purityEnrichment.processVariant(x.context()));
    }

    private boolean isFittingCandidate(final SomaticVariant variant)
    {
        if(variant.type() != VariantType.SNP)
            return false;

        if(!variant.isPass())
            return false;

        if(!mFilter.test(variant.context()))
            return false;

        if(!hasTumorDepth(variant))
            return false;

        return true;
    }

    private boolean hasTumorDepth(final SomaticVariant variant)
    {
        if(!mConfig.runTumor() || !variant.hasAlleleDepth())
            return false;

        return variant.tumorAlleleDepth().totalReadCount() > 0;
    }
}
