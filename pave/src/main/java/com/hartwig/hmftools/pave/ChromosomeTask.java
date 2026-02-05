package com.hartwig.hmftools.pave;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.pon.GnomadCache.PON_GNOMAD_FILTER;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILTER;
import static com.hartwig.hmftools.pave.FilterType.PANEL;
import static com.hartwig.hmftools.pave.FilterType.PASS;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.pave.PaveConfig.isSbx;
import static com.hartwig.hmftools.pave.PaveConfig.isUltima;
import static com.hartwig.hmftools.pave.PaveConstants.GNMOAD_FILTER_HOTSPOT_PATHOGENIC_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.GNMOAD_FILTER_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_REPEAT_COUNT;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_SBX_FACTOR;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_SBX_MIN_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_ULTIMA_FACTOR;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_ULTIMA_MIN_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_INDEL_ARTEFACT_VAF_REDUCTION;
import static com.hartwig.hmftools.pave.PaveConstants.PON_MEAN_READ_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_REPEAT_COUNT_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_SAMPLE_COUNT_THRESHOLD;
import static com.hartwig.hmftools.pave.PaveConstants.PON_VAF_THRESHOLD;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.annotation.PonAnnotation.PON_ARTEFACT_FILTER;
import static com.hartwig.hmftools.pave.impact.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.impact.PaveUtils.findVariantImpacts;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.pon.GnomadChrCache;
import com.hartwig.hmftools.common.variant.pon.MultiPonStatus;
import com.hartwig.hmftools.common.variant.pon.PonChrCache;
import com.hartwig.hmftools.common.variant.pon.PonVariantData;
import com.hartwig.hmftools.pave.annotation.ClinvarChrCache;
import com.hartwig.hmftools.pave.annotation.GnomadAnnotation;
import com.hartwig.hmftools.pave.annotation.MappabilityChrCache;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;
import com.hartwig.hmftools.pave.annotation.ReferenceData;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantImpactBuilder;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomeTask implements Callable<Void>
{
    private final HumanChromosome mChromosome;
    private final String mChromosomeStr;
    private final PaveConfig mConfig;
    private final VcfWriter mVcfWriter;
    private final TranscriptWriter mTranscriptWriter;

    private final ReferenceData mReferenceData;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;

    // local chromosome annotation caches
    private GnomadChrCache mGnomadCache;
    private ClinvarChrCache mClinvarCache;
    private MappabilityChrCache mMappability;
    private PonChrCache mStandardPon;
    private PonChrCache mArtefactsPon;
    private final GeneCacheIndexing mGeneCacheIndexing;

    public ChromosomeTask(
            final HumanChromosome chromosome, final PaveConfig config, final ReferenceData referenceData,
            final VcfWriter vcfWriter, final TranscriptWriter transcriptWriter)
    {
        mChromosome = chromosome;
        mChromosomeStr = config.RefGenVersion.versionedChromosome(chromosome.toString());
        mConfig = config;
        mVcfWriter = vcfWriter;
        mTranscriptWriter = transcriptWriter;
        mReferenceData = referenceData;

        mImpactBuilder = new VariantImpactBuilder(mReferenceData.GeneDataCache);
        mImpactClassifier = new ImpactClassifier(mReferenceData.RefGenome);

        mGnomadCache = null;
        mClinvarCache = null;
        mMappability = null;
        mStandardPon = null;
        mArtefactsPon = null;
        mGeneCacheIndexing = mReferenceData.GeneDataCache.createIndexing(mChromosomeStr);
    }

    @Override
    public Void call()
    {
        int variantCount = 0;

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.VcfFile, mConfig.requireIndex());

        if(!vcfFileReader.fileValid())
        {
            PV_LOGGER.error("invalid somatic VCF file({})", mConfig.VcfFile);
            System.exit(1);
        }

        mGnomadCache = mReferenceData.Gnomad.getChromosomeCache(mChromosomeStr);
        mClinvarCache = mReferenceData.Clinvar.getChromosomeCache(mChromosomeStr);
        mMappability = mReferenceData.VariantMappability.getChromosomeCache(mChromosomeStr);
        mStandardPon = mReferenceData.StandardPon.getChromosomeCache(mChromosomeStr);
        mArtefactsPon = mReferenceData.ArtefactsPon.getChromosomeCache(mChromosomeStr);

        ChrBaseRegion chrRegion;

        if(mConfig.SpecificChrRegions.Regions.size() == 1 && mConfig.SpecificChrRegions.Regions.get(0).Chromosome.equals(mChromosomeStr))
        {
            chrRegion = mConfig.SpecificChrRegions.Regions.get(0);
        }
        else
        {
            RefGenomeCoordinates coordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
            chrRegion = new ChrBaseRegion(mChromosomeStr, 1, coordinates.Lengths.get(mChromosome));
        }

        PV_LOGGER.debug("chr({}) starting variant annotation", mChromosomeStr);

        CloseableTribbleIterator<VariantContext> varIterator = mConfig.requireIndex() ?
                vcfFileReader.regionIterator(chrRegion) : vcfFileReader.iterator();

        for(VariantContext variantContext : varIterator)
        {
            if(mConfig.SpecificChrRegions.hasFilters())
            {
                if(!mConfig.SpecificChrRegions.Regions.isEmpty()
                && mConfig.SpecificChrRegions.Regions.stream().noneMatch(x -> x.containsPosition(
                        variantContext.getContig(), variantContext.getStart())))
                {
                    continue;
                }
            }

            processVariant(variantContext);
            ++variantCount;

            if(variantCount > 0 && (variantCount % 100000) == 0)
            {
                PV_LOGGER.debug("chr({}) processed {} variants", mChromosome, variantCount);
            }
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET);

        PV_LOGGER.info("chr({}) complete for {} variants", mChromosome, variantCount);

        mReferenceData.onChromosomeComplete(mChromosomeStr);

        mVcfWriter.onChromosomeComplete(mChromosome);

        return null;
    }

    private void processVariant(final VariantContext variantContext)
    {
        if(!HumanChromosome.contains(variantContext.getContig()))
            return;

        VariantData variant = VariantData.fromContext(variantContext);

        boolean isPass = variantContext.getFilters().isEmpty() || variantContext.getFilters().contains(PASS);

        if(mConfig.Filter == FilterType.PASS)
        {
            if(!isPass)
                return;
        }
        else if(mConfig.Filter == PANEL)
        {
            VariantTier tier = VariantTier.fromContext(variantContext);

            // anything in the panel or passing variants
            if(tier != VariantTier.HOTSPOT && tier != VariantTier.PANEL && !isPass)
                return;
        }

        try
        {
            variant.setRealignedVariant(createRightAlignedVariant(variant, mImpactClassifier.refGenome()));

            findVariantImpacts(variant, mImpactClassifier, mReferenceData.GeneDataCache, mGeneCacheIndexing);

            processPhasedVariants(variant.localPhaseSet());

            if(!variant.hasLocalPhaseSet())
                processVariant(variant);
        }
        catch(Exception e)
        {
            PV_LOGGER.error("error processing var({})", variant);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processPhasedVariants(int currentLocalPhaseSet)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
            variants.forEach(x -> processVariant(x));
    }

    private void processVariant(final VariantData variant)
    {
        // can be null if no impacts exist for any transcript
        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        annotateAndFilter(variant);

        if(mConfig.SetReportable)
            mReferenceData.ReportableClassifier.setReportability(variant, variantImpact);

        if(mConfig.WritePassOnly && !variant.filters().isEmpty())
            return;

        VariantContext newVariant = mVcfWriter.buildVariant(variant.context(), variant, variantImpact);
        mVcfWriter.writeVariant(mChromosome, newVariant);

        if(mConfig.WriteTranscriptFile)
        {
            for(Map.Entry<String, List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
            {
                final String geneName = entry.getKey();
                mTranscriptWriter.writeVariantData(variant, geneName);
            }
        }
    }

    private void annotateAndFilter(final VariantData variant)
    {
        if(mMappability != null)
            mMappability.annotateVariant(variant);

        if(mClinvarCache != null)
            mClinvarCache.annotateVariant(variant);

        mReferenceData.BlacklistedVariants.annotateVariant(variant);

        if(mGnomadCache != null)
            mReferenceData.Gnomad.annotateVariant(variant, mGnomadCache);
        else
            GnomadAnnotation.annotateFromContext(variant);

        if(mStandardPon != null)
            mReferenceData.StandardPon.annotateVariant(variant, mStandardPon);
        else
            PonAnnotation.annotateFromContext(variant);

        applyFilters(variant);
    }

    private void applyFilters(final VariantData variant)
    {
        applyFilters(variant, mConfig.SampleId, mReferenceData.StandardPon, mArtefactsPon, mReferenceData.Gnomad.applyFilter());
    }

    @VisibleForTesting
    public static void applyFilters(
            final VariantData variant, final String sampleId, final PonAnnotation standardPon, final PonChrCache artefactsPon,
            boolean applyGnomadFilter)
    {
        variant.filters().clear();

        boolean isHotspot = variant.tier() == VariantTier.HOTSPOT;
        boolean clinvarPathogenic = PathogenicSummaryFactory.fromContext(variant.context()).Status.isPathogenic();
        boolean hotspotOrPathogenic = isHotspot || clinvarPathogenic;
        int repeatCount = variant.repeatCount();

        // the WGS PON
        boolean ponFilter = standardPon.filterOnTierCriteria(variant.tier(), variant.ponSampleCount(), variant.ponMaxReadCount());

        if(hotspotOrPathogenic)
        {
            if(variant.ponSampleCount() == 0 || belowPonThreshold(variant.ponSampleCount(), repeatCount))
            {
                ponFilter = false;
            }
            else
            {
                double variantVaf = variant.sampleVaf(sampleId);
                int repeatBaseLength = variant.repeatSequence().length() * repeatCount;
                double vafLimit = max(PON_VAF_THRESHOLD, 0.01 * repeatBaseLength);

                if(variant.ponMeanReadCount() < PON_MEAN_READ_THRESHOLD && variantVaf > vafLimit)
                {
                    ponFilter = false;
                }
            }
        }

        if(ponFilter && variant.isIndel() && repeatCount >= PON_INDEL_ARTEFACT_REPEAT_COUNT
        && (isSbx() || isUltima()) && variant.multiPonStatus() == MultiPonStatus.ARTEFACT)
        {
            // the AF-based indel filtering would only be applied to variants with SeqTechOnly=true
            //	- Variant is only in seq tech PON (not in Illumina PON)
            //	- Variant is MSI indel with repeat count >= 7
            //	- Variant has PON_MAX < min(60 * (sampleAF - 0.1), 18) (if Ultima) or min(40 * (sampleAF - 0.1), 8) (if SBX)
            double variantVaf = variant.sampleVaf(sampleId);
            int ponMinThreshold, ponFactor;

            if(isUltima())
            {
                ponMinThreshold = PON_INDEL_ARTEFACT_ULTIMA_MIN_THRESHOLD;
                ponFactor = PON_INDEL_ARTEFACT_ULTIMA_FACTOR;
            }
            else
            {
                ponMinThreshold = PON_INDEL_ARTEFACT_SBX_MIN_THRESHOLD;
                ponFactor = PON_INDEL_ARTEFACT_SBX_FACTOR;
            }

            double ponThreshold = min(ponFactor * (variantVaf - PON_INDEL_ARTEFACT_VAF_REDUCTION), ponMinThreshold);

            if(variant.ponSampleCount() < ponThreshold)
                ponFilter = false;
        }

        if(ponFilter)
            variant.addFilter(PON_FILTER);

        if(artefactsPon != null)
        {
            PonVariantData artefactPonData = artefactsPon.getPonData(variant.Position, variant.Ref, variant.Alt);

            if(artefactPonData != null)
            {
                if(!hotspotOrPathogenic || !belowPonThreshold(artefactPonData.Samples, repeatCount))
                {
                    variant.addFilter(PON_ARTEFACT_FILTER);
                }
            }
        }

        if(applyGnomadFilter)
        {
            Double gnmoadFrequency = variant.gnomadFrequency();
            double gnomadThreshold = hotspotOrPathogenic ? GNMOAD_FILTER_HOTSPOT_PATHOGENIC_THRESHOLD : GNMOAD_FILTER_THRESHOLD;

            if(gnmoadFrequency != null && gnmoadFrequency >= gnomadThreshold)
            {
                variant.addFilter(PON_GNOMAD_FILTER);
            }
        }
    }

    private static boolean belowPonThreshold(final int ponSampleCount, final int repeatCount)
    {
        return ponSampleCount < PON_SAMPLE_COUNT_THRESHOLD && repeatCount < PON_REPEAT_COUNT_THRESHOLD;
    }
}
