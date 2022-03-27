package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact.VAR_TRANS_IMPACT_ANNOATATION;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.purple.SomaticVariantCache;
import com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
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
        mOutputVCF = config.OutputDir + config.TumorId + PURPLE_SOMATIC_VCF_SUFFIX;
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

        mVcfWriter = null;
    }

    public double microsatelliteIndelsPerMb()
    {
        return mMicrosatelliteIndels.microsatelliteIndelsPerMb();
    }

    public MicrosatelliteStatus microsatelliteStatus()
    {
        return mEnabled ? MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb()) : MicrosatelliteStatus.UNKNOWN;
    }

    public double tumorMutationalBurdenPerMb()
    {
        return mTumorMutationalLoad.burdenPerMb();
    }

    public int tumorMutationalLoad()
    {
        return mTumorMutationalLoad.load();
    }

    public TumorMutationalStatus tumorMutationalBurdenPerMbStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromBurdenPerMb(tumorMutationalBurdenPerMb()) : TumorMutationalStatus.UNKNOWN;
    }

    public TumorMutationalStatus tumorMutationalLoadStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromLoad(tumorMutationalLoad()) : TumorMutationalStatus.UNKNOWN;
    }

    public List<DriverCatalog> buildDrivers(final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        return mDrivers.buildCatalog(geneCopyNumberMap);
    }

    public Set<String> reportedGenes()
    {
        return mReportedGenes;
    }

    public List<VariantContextDecorator> downsampledVariants() { return mDownsampledVariants; }

    public void processAndWrite(
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions)
    {
        if(!mEnabled || mPeakModel == null)
            return;

        try
        {
            VCFHeader readHeader = mSomaticVariants.getVcfHeader();
            boolean isPaveAnnotated = readHeader.hasInfoLine(VAR_TRANS_IMPACT_ANNOATATION);

            if(!isPaveAnnotated)
            {
                PPL_LOGGER.info("SnpEff annotation enabled");
            }

            mVcfWriter = new VariantContextWriterBuilder().setOutputFile(mOutputVCF)
                    .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final SomaticVariantEnrichment enricher = new SomaticVariantEnrichment(
                    !isPaveAnnotated, mConfig.SomaticFitting.clonalityBinWidth(), mConfig.Version,
                    mConfig.ReferenceId, mConfig.TumorId, mReferenceData, purityAdjuster, copyNumbers, fittedRegions, mPeakModel);

            final VCFHeader header = enricher.populateHeader(readHeader);
            mVcfWriter.writeHeader(header);

            boolean tumorOnly = mConfig.tumorOnlyMode();

            int flushCount = 10000;
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
                }
            }

            enricher.flush(); // finalise any enrichment routines with queued variants

            // write enriched variants to VCF
            for(SomaticVariant variant : mSomaticVariants.variants())
            {
                boolean isValidChromosome = HumanChromosome.contains(variant.chromosome());

                if(isValidChromosome && variant.isPass())
                {
                    mTumorMutationalLoad.processVariant(variant);
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
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to enrich somatic variants: {}", e.toString());
        }
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
