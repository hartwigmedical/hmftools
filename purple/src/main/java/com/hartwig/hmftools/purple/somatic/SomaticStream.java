package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.CodingEffect.hasProteinImpact;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_MISSENSE;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.ASSUMED_BIALLELIC_FRACTION;
import static com.hartwig.hmftools.purple.PurpleConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.addReportableTranscriptList;
import static com.hartwig.hmftools.purple.somatic.SomaticVariantEnrichment.populateHeader;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.purple.fittingsnv.PeakModelData;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.plot.RChartData;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticStream
{
    private final ReferenceData mReferenceData;
    private final PurpleConfig mConfig;

    private boolean mEnabled;
    private final String mOutputVCF;
    private final TumorMutationalLoad mTumorMutationalLoad;
    private final SomaticGermlineLikelihood mSomaticGermlineLikelihood;
    private final MicrosatelliteIndels mMicrosatelliteIndels;
    private final SomaticVariantDrivers mDrivers;
    private final SomaticVariantCache mSomaticVariants;
    private final RChartData mRChartData;
    private final DriverGenePanel mGenePanel;
    private final List<PeakModelData> mPeakModelData;
    private final Set<String> mReportedGenes;

    private final List<VariantContextDecorator> mDownsampledVariants; // cached for charting
    private final int mSnpMod;
    private final int mIndelMod;
    private int mSnpCount;
    private int mIndelCount;

    // computed values
    private double mTml;
    private double mTmb;
    private double mMsiIndelPerMb;

    private VariantContextWriter mVcfWriter;

    private static final int CHART_DOWNSAMPLE_FACTOR = 25000; // eg for 50K variants, only every second will be kept for plotting

    public SomaticStream(
            final PurpleConfig config, final ReferenceData referenceData, final SomaticVariantCache somaticVariantCache)
    {
        mReferenceData = referenceData;
        mConfig = config;

        mGenePanel = referenceData.DriverGenes;
        mSomaticVariants = somaticVariantCache;
        mPeakModelData = Lists.newArrayList();
        mOutputVCF = PurpleCommon.purpleSomaticVcfFile(config.OutputDir, config.TumorId);
        mEnabled = somaticVariantCache.hasData();
        mTumorMutationalLoad = new TumorMutationalLoad(mReferenceData.TargetRegions);
        mSomaticGermlineLikelihood = new SomaticGermlineLikelihood(mConfig, somaticVariantCache.genotypeIds());
        mMicrosatelliteIndels = new MicrosatelliteIndels(mReferenceData.TargetRegions);
        mDrivers = new SomaticVariantDrivers(mGenePanel);
        mRChartData = new RChartData(config, config.TumorId);

        mReportedGenes = Sets.newHashSet();
        mDownsampledVariants = Lists.newArrayList();
        mSnpMod = somaticVariantCache.snpCount() <= CHART_DOWNSAMPLE_FACTOR ? 1 : somaticVariantCache.snpCount() / CHART_DOWNSAMPLE_FACTOR;
        mIndelMod = somaticVariantCache.indelCount() <= CHART_DOWNSAMPLE_FACTOR ? 1 : somaticVariantCache.indelCount() / CHART_DOWNSAMPLE_FACTOR;

        mTmb = 0;
        mTml = 0;
        mMsiIndelPerMb = 0;

        mVcfWriter = null;
    }

    public double msiIndelsPerMb() { return mMsiIndelPerMb; }
    public double tumorMutationalBurdenPerMb() { return mTmb; }
    public int tumorMutationalLoad() { return (int)round(mTml); }

    public MicrosatelliteStatus microsatelliteStatus()
    {
        return mEnabled ? MicrosatelliteStatus.fromIndelsPerMb(mMsiIndelPerMb) : MicrosatelliteStatus.UNKNOWN;
    }

    public TumorMutationalStatus tumorMutationalBurdenPerMbStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromBurdenPerMb(mTmb) : TumorMutationalStatus.UNKNOWN;
    }

    public TumorMutationalStatus tumorMutationalLoadStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromLoad(mTml) : TumorMutationalStatus.UNKNOWN;
    }

    public List<DriverCatalog> buildDrivers(final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        if(mReferenceData.TargetRegions.hasTargetRegions())
        {
            // override the counts passed to the DNDS calcs
            int inferredSnvCount = (int)round(mTml * mReferenceData.TargetRegions.tmbRatio() * MB_PER_GENOME);
            int inferredIndelCount = (int)round(mMsiIndelPerMb * MB_PER_GENOME);

            Map<VariantType,Integer> variantTypeCounts = Maps.newHashMap();
            variantTypeCounts.put(VariantType.SNP, inferredSnvCount);
            variantTypeCounts.put(VariantType.INDEL, inferredIndelCount);

            Map<VariantType,Integer> variantTypeCountsbiallelic = Maps.newHashMap();
            variantTypeCountsbiallelic.put(VariantType.SNP, (int)(inferredSnvCount * ASSUMED_BIALLELIC_FRACTION));
            variantTypeCountsbiallelic.put(VariantType.INDEL, (int)(inferredIndelCount * ASSUMED_BIALLELIC_FRACTION));

            PPL_LOGGER.debug("target-regions inferred driver variant counts snv({}) indel({})",
                    inferredSnvCount, inferredIndelCount);

            mDrivers.overrideVariantCounts(variantTypeCounts, variantTypeCountsbiallelic);
        }

        return mDrivers.buildCatalog(geneCopyNumberMap);
    }

    public Set<String> reportedGenes() { return mReportedGenes; }
    public List<VariantContextDecorator> downsampledVariants() { return mDownsampledVariants; }
    public List<PeakModelData> peakModelData() { return mPeakModelData; }

    public void processAndWrite(final PurityAdjuster purityAdjuster)
    {
        if(!mEnabled)
            return;

        PPL_LOGGER.debug("modelling somatic peaks");
        SomaticPeakStream somaticPeakStream = new SomaticPeakStream();

        mPeakModelData.addAll(somaticPeakStream.generateModelPeaks(mSomaticVariants));

        try
        {
            VCFHeader vcfHeader = mSomaticVariants.getVcfHeader();

            mVcfWriter = new VariantContextWriterBuilder().setOutputFile(mOutputVCF)
                    .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            populateHeader(vcfHeader, mConfig.Version);

            if(mConfig.tumorOnlyMode())
                SomaticGermlineLikelihood.enrichHeader(vcfHeader);

            mVcfWriter.writeHeader(vcfHeader);

            boolean tumorOnly = mConfig.tumorOnlyMode();
            AtomicInteger kataegisId = new AtomicInteger();

            if(mConfig.Threads > 1)
            {
                List<SomaticVariantEnrichment> enrichers = Lists.newArrayList();

                for(int i = 0; i < mConfig.Threads; ++i)
                {
                    enrichers.add(new SomaticVariantEnrichment(i, mConfig, mReferenceData, mPeakModelData, kataegisId));
                }

                int taskIndex = 0;
                String currentChr = !mSomaticVariants.variants().isEmpty() ? mSomaticVariants.variants().get(0).chromosome() : "";

                for(SomaticVariant variant : mSomaticVariants.variants())
                {
                    if(!currentChr.equals(variant.chromosome()))
                    {
                        currentChr = variant.chromosome();
                        ++taskIndex;

                        if(taskIndex >= enrichers.size())
                            taskIndex = 0;
                    }

                    enrichers.get(taskIndex).addVariant(variant);
                }

                final List<Callable> callableList = enrichers.stream().collect(Collectors.toList());
                TaskExecutor.executeTasks(callableList, mConfig.Threads);
            }
            else
            {
                SomaticVariantEnrichment enricher = new SomaticVariantEnrichment(0, mConfig, mReferenceData, mPeakModelData, kataegisId);
                mSomaticVariants.variants().forEach(x -> enricher.addVariant(x));
                enricher.call();
            }

            // various processing for charting, TMB/L calcs, drivers
            for(SomaticVariant variant : mSomaticVariants.variants())
            {
                if(!HumanChromosome.contains(variant.chromosome()))
                    continue;

                if(variant.isPass() || mConfig.WriteAllSomatics)
                    mSomaticGermlineLikelihood.processVariant(variant, purityAdjuster.purity());

                if(variant.isPass())
                {
                    mTumorMutationalLoad.processVariant(variant);
                    mMicrosatelliteIndels.processVariant(variant);
                    checkDrivers(variant, true); // sets reportable flag if applicable

                    mRChartData.processVariant(variant);
                    checkChartDownsampling(variant);
                }
            }

            // should not be required if coding effects have been set correctly for phased variants in Pave
            checkPhasedReportableVariants();

            // write enriched variants to VCF
            for(SomaticVariant variant : mSomaticVariants.variants())
            {
                if(!tumorOnly || variant.isPass() || mConfig.WriteAllSomatics)
                    mVcfWriter.add(variant.context());
            }

            mVcfWriter.close();
            mRChartData.write();

            calculateVariantLoadValues();

            PPL_LOGGER.debug("charting variants: total(snvs={} indels={}) downsampled({} snvMod={} indelMod={})",
                    mSnpCount, mIndelCount, mDownsampledVariants.size(), mSnpMod, mIndelMod);
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to enrich somatic variants: {}", e.toString());
        }
    }

    private void calculateVariantLoadValues()
    {
        mTml = mTumorMutationalLoad.calcTml();
        mMsiIndelPerMb = mMicrosatelliteIndels.calcMsiIndelsPerMb();

        if(mReferenceData.TargetRegions.hasTargetRegions())
        {
            mTmb = mMsiIndelPerMb + mTml * mReferenceData.TargetRegions.tmbRatio();
        }
        else
        {
            mTmb = mTumorMutationalLoad.burden() / MB_PER_GENOME;
        }

        PPL_LOGGER.info(String.format("load(%.1f tml=%.4f) msiIndels(%d perMb=%.4f) burden(%.1f perMb=%.4f)",
                mTumorMutationalLoad.load(), mTml, mMicrosatelliteIndels.msiIndelCount(), mMsiIndelPerMb,
                mTumorMutationalLoad.burden(), mTmb));
    }

    private void checkDrivers(final SomaticVariant variant, boolean updateVcf)
    {
        boolean reported = mDrivers.checkSomaticVariant(variant);

        if(reported && updateVcf)
        {
            variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

            mReportedGenes.add(variant.decorator().gene());

            // check alt transcript status vs canonical
            addReportableTranscriptList(variant.type(), variant.context(), variant.variantImpact());
        }
    }

    public void registerReportedVariants()
    {
        for(SomaticVariant variant : mSomaticVariants.variants())
        {
            boolean isValidChromosome = HumanChromosome.contains(variant.chromosome());

            if(isValidChromosome && variant.isPass())
            {
                checkDrivers(variant, false);
            }
        }
    }

    private static boolean hasPhasedEffect(final SomaticVariant variant)
    {
        return variant.variantImpact().CanonicalEffect.contains(PHASED_INFRAME_INSERTION.effect())
            || variant.variantImpact().CanonicalEffect.contains(PHASED_INFRAME_DELETION.effect())
            || variant.variantImpact().CanonicalEffect.contains(PHASED_MISSENSE.effect());
    }

    private void checkPhasedReportableVariants()
    {
        // any non-reportable variant that forms a phased inframe INDEL with a reportable variant is marked as reportable too
        for(int i = 0; i < mSomaticVariants.variants().size(); ++i)
        {
            SomaticVariant variant = mSomaticVariants.variants().get(i);

            // first find any reportable phased inframe INDEL
            if(!variant.context().hasAttribute(REPORTED_FLAG) || !hasPhasedEffect(variant))
                continue;

            List<Integer> localPhaseSets = variant.context().getAttributeAsIntList(LOCAL_PHASE_SET, 0);

            if(localPhaseSets.isEmpty())
                continue;

            // look forwards and backwards for unreported passing variants in the same phase set
            for(int direction = 0; direction <= 1; ++direction)
            {
                boolean searchBack = direction == 0;

                int j = i;
                while(true)
                {
                    if(searchBack)
                        --j;
                    else
                        ++j;

                    if(j < 0 || j >= mSomaticVariants.variants().size())
                        break;

                    SomaticVariant nextVariant = mSomaticVariants.variants().get(j);

                    if(!nextVariant.isPass() || nextVariant.context().hasAttribute(REPORTED_FLAG) || !hasPhasedEffect(variant))
                        continue;

                    // must have a coding impact
                    if(nextVariant.variantImpact() == null || !hasProteinImpact(nextVariant.variantImpact().CanonicalCodingEffect))
                        continue;

                    List<Integer> nextLocalPhaseSets = nextVariant.context().getAttributeAsIntList(LOCAL_PHASE_SET, 0);

                    // stop looking when phase set changes or is empty, so assumes that there aren't unphased variants in between
                    if(nextLocalPhaseSets.isEmpty())
                        break;

                    if(nextLocalPhaseSets.stream().noneMatch(x -> localPhaseSets.contains(x)))
                        break;

                    nextVariant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

                    PPL_LOGGER.debug("var({}) setting reported due to inframe-phasing with other({})", nextVariant, variant);

                    // add to appropriate driver caches for DNDS calcs
                    mDrivers.addPhasedReportableVariant(nextVariant, variant);
                }
            }
        }
    }

    private void checkChartDownsampling(final SomaticVariant variant)
    {
        if(mConfig.Charting.Disabled)
            return;

        if(!HumanChromosome.contains(variant.chromosome()))
            return;

        if(variant.type() == VariantType.INDEL)
        {
            mIndelCount++;

            if(mIndelCount % mIndelMod == 0)
                mDownsampledVariants.add(variant.decorator());
        }
        else
        {
            mSnpCount++;

            if(mSnpCount % mSnpMod == 0)
                mDownsampledVariants.add(variant.decorator());
        }
    }
}
