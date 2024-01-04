package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.impact.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.impact.PaveUtils.findVariantImpacts;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.VcfWriter.buildVariant;
import static com.hartwig.hmftools.pave.annotation.PonAnnotation.PON_ARTEFACT_FILTER;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.annotation.ClinvarChrCache;
import com.hartwig.hmftools.pave.annotation.GnomadChrCache;
import com.hartwig.hmftools.pave.annotation.MappabilityChrCache;
import com.hartwig.hmftools.pave.annotation.PonChrCache;
import com.hartwig.hmftools.pave.annotation.ReferenceData;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantImpactBuilder;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomeTask implements Callable
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
    public Long call()
    {
        int variantCount = 0;

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.VcfFile, true);

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

        RefGenomeCoordinates coordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        ChrBaseRegion chrRegion = new ChrBaseRegion(mChromosomeStr, 1, coordinates.Lengths.get(mChromosome));

        PV_LOGGER.debug("chr({}) starting variant annotation", mChromosome);

        for(VariantContext variantContext : vcfFileReader.regionIterator(chrRegion))
        {
            if(!mConfig.SpecificRegions.isEmpty())
            {
                if(mConfig.SpecificRegions.stream().noneMatch(x -> x.containsPosition(variantContext.getContig(), variantContext.getStart())))
                    continue;
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

        return (long)0;
    }

    private void processVariant(final VariantContext variantContext)
    {
        if(!HumanChromosome.contains(variantContext.getContig()))
            return;

        VariantData variant = VariantData.fromContext(variantContext);

        if(mConfig.ReadPassOnly)
        {
            if(!variantContext.getFilters().isEmpty() && !variantContext.getFilters().contains(PASS_FILTER))
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

        VariantContext newVariant = buildVariant(variant.context(), variant, variantImpact);
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
        if(mGnomadCache != null)
            mReferenceData.Gnomad.annotateVariant(variant, mGnomadCache);

        if(mMappability != null)
            mMappability.annotateVariant(variant);

        if(mClinvarCache != null)
            mClinvarCache.annotateVariant(variant);

        mReferenceData.BlacklistedVariants.annotateVariant(variant);

        if(mStandardPon != null)
            mReferenceData.StandardPon.annotateVariant(variant, mStandardPon);

        if(mArtefactsPon != null && mArtefactsPon.getPonData(variant) != null)
            variant.addFilter(PON_ARTEFACT_FILTER);
    }
}
