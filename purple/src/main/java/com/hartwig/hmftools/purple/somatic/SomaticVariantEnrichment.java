package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CLONALITY_BIN_WIDTH;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.PurpleVcfTags;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.fitting.PeakModelData;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;

import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment implements Callable
{
    private final KataegisEnrichment mKataegisEnrichment;
    private final SomaticRefContextEnrichment mSomaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment mSubclonalLikelihoodEnrichment;
    private final SomaticGenotypeEnrichment mGenotypeEnrichment;

    private final PurpleConfig mConfig;
    private final int mTaskId;
    private final List<SomaticVariant> mVariants;

    public SomaticVariantEnrichment(
            final int taskId, final PurpleConfig config, final ReferenceData refData, final List<PeakModelData> peakModel,
            final AtomicInteger kataegisId)
    {
        mConfig = config;
        mTaskId = taskId;
        mGenotypeEnrichment = new SomaticGenotypeEnrichment(mConfig.ReferenceId, mConfig.TumorId);
        mSubclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(CLONALITY_BIN_WIDTH, peakModel);
        mKataegisEnrichment = new KataegisEnrichment(kataegisId);
        mSomaticRefContextEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, null);
        mVariants = Lists.newArrayList();
    }

    public void addVariant(final SomaticVariant variant) { mVariants.add(variant); }

    @Override
    public Long call()
    {
        int flushCount = 100000;
        // int gcCount = 250000;
        int varCount = 0;

        boolean tumorOnly = mConfig.tumorOnlyMode();

        for(SomaticVariant variant : mVariants)
        {
            if(tumorOnly && variant.isFiltered() && !mConfig.WriteAllSomatics)
                continue;

            enrich(variant);
            ++varCount;

            if(varCount > 0 && (varCount % flushCount) == 0)
            {
                PPL_LOGGER.debug("{}: enriched {} somatic variants", mTaskId, varCount);
            }
        }

        flush(); // finalise any enrichment routines with queued variants

        return (long)0;
    }

    public void enrich(final SomaticVariant variant)
    {
        mSomaticRefContextEnrichment.processVariant(variant.context());

        mKataegisEnrichment.processVariant(variant);

        mSubclonalLikelihoodEnrichment.processVariant(variant);

        mGenotypeEnrichment.processVariant(variant);
    }

    public void flush()
    {
        mKataegisEnrichment.flush();
    }

    public static VCFHeader populateHeader(final VCFHeader template, final String purpleVersion)
    {
        VCFHeader header = SageVcfTags.addRefContextHeader(template);
        header = KataegisEnrichment.enrichHeader(header);
        header = SubclonalLikelihoodEnrichment.enrichHeader(header);
        header = HotspotEnrichment.enrichHeader(header);
        header = PurpleVcfTags.addSomaticHeader(purpleVersion, header);
        return header;
    }
}
