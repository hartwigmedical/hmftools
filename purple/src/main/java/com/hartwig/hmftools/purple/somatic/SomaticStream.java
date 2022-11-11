package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MB_PER_GENOME;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.purple.SomaticVariantCache;
import com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.purple.config.PurpleConfig;
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
    private final MicrosatelliteIndels mMicrosatelliteIndels;
    private final SomaticVariantDrivers mDrivers;
    private final SomaticVariantCache mSomaticVariants;
    private final RChartData mRChartData;
    private final DriverGenePanel mGenePanel;
    private final List<PeakModel> mPeakModel;
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
            final PurpleConfig config, final ReferenceData referenceData, final SomaticVariantCache somaticVariants,
            final List<PeakModel> peakModel)
    {
        mReferenceData = referenceData;
        mConfig = config;

        mGenePanel = referenceData.DriverGenes;
        mPeakModel = peakModel;
        mOutputVCF = PurpleCommon.purpleSomaticVcfFile(config.OutputDir, config.TumorId);
        mEnabled = somaticVariants.hasData();
        mTumorMutationalLoad = new TumorMutationalLoad(mReferenceData.TargetRegions);
        mMicrosatelliteIndels = new MicrosatelliteIndels(mReferenceData.TargetRegions);
        mDrivers = new SomaticVariantDrivers(mGenePanel);
        mSomaticVariants = somaticVariants;
        mRChartData = new RChartData(config.OutputDir, config.TumorId);

        mReportedGenes = Sets.newHashSet();
        mDownsampledVariants = Lists.newArrayList();
        mSnpMod = somaticVariants.snpCount() <= CHART_DOWNSAMPLE_FACTOR ? 1 : somaticVariants.snpCount() / CHART_DOWNSAMPLE_FACTOR;
        mIndelMod = somaticVariants.indelCount() <= CHART_DOWNSAMPLE_FACTOR ? 1 : somaticVariants.indelCount() / CHART_DOWNSAMPLE_FACTOR;

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
        return mDrivers.buildCatalog(geneCopyNumberMap);
    }

    public Set<String> reportedGenes() { return mReportedGenes; }

    public List<VariantContextDecorator> downsampledVariants() { return mDownsampledVariants; }

    public void processAndWrite(
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
        if(!mEnabled || mPeakModel == null)
            return;

        try
        {
            VCFHeader readHeader = mSomaticVariants.getVcfHeader();

            mVcfWriter = new VariantContextWriterBuilder().setOutputFile(mOutputVCF)
                    .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final SomaticVariantEnrichment enricher = new SomaticVariantEnrichment(
                    mConfig.Version, mConfig.ReferenceId, mConfig.TumorId, mReferenceData, purityAdjuster, copyNumbers, fittedRegions, mPeakModel);

            final VCFHeader header = enricher.populateHeader(readHeader);
            mVcfWriter.writeHeader(header);

            boolean tumorOnly = mConfig.tumorOnlyMode();

            int flushCount = 100000;
            int gcCount = 250000;
            int varCount = 0;

            for(SomaticVariant variant : mSomaticVariants.variants())
            {
                if(tumorOnly && variant.isFiltered())
                    continue;

                enricher.enrich(variant);
                ++varCount;

                if(varCount > 0 && (varCount % flushCount) == 0)
                {
                    PPL_LOGGER.debug("enriched {} somatic variants", varCount);

                    if((varCount % gcCount) == 0)
                        System.gc();
                }
            }

            enricher.flush(); // finalise any enrichment routines with queued variants

            // write enriched variants to VCF
            for(SomaticVariant variant : mSomaticVariants.variants())
            {
                boolean isValidChromosome = HumanChromosome.contains(variant.chromosome());

                if(isValidChromosome && variant.isPass())
                {
                    mTumorMutationalLoad.processVariant(variant, purityAdjuster.purity());
                    mMicrosatelliteIndels.processVariant(variant);
                    checkDrivers(variant, true); // sets reportable flag if applicable

                    mRChartData.processVariant(variant);
                    checkChartDownsampling(variant);
                }

                // expect only pass or PON to be loaded into Purple - in tumor-only mode, the PON variants should be dropped
                if(!tumorOnly || variant.isPass())
                    mVcfWriter.add(variant.context());
            }

            mVcfWriter.close();
            mRChartData.write();

            calculateVariantLoadValues();
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
        }
    }

    private void checkChartDownsampling(final SomaticVariant variant)
    {
        if(mConfig.Charting.disabled())
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
