package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.CLONALITY_BIN_WIDTH;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.PurpleVcfTags;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.fittingsnv.PeakModelData;
import com.hartwig.hmftools.purple.ReferenceData;

import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment implements Callable
{
    private final KataegisEnrichment mKataegisEnrichment;
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
        mKataegisEnrichment.processVariant(variant);
        mSubclonalLikelihoodEnrichment.processVariant(variant);
        mGenotypeEnrichment.processVariant(variant);
    }

    public void flush()
    {
        mKataegisEnrichment.flush();
    }

    public static void populateHeader(final VCFHeader header, final String purpleVersion)
    {
        KataegisEnrichment.enrichHeader(header);
        SubclonalLikelihoodEnrichment.enrichHeader(header);
        HotspotEnrichment.enrichHeader(header);
        PurpleVcfTags.addSomaticHeader(purpleVersion, header);
    }
}
