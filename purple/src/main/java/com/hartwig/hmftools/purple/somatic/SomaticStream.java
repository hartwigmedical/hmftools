package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.SomaticVariantDrivers;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalLoad;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.plot.RChartData;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticStream
{
    private static final int MAX_DOWNSAMPLE = 25000;

    private final ReferenceData mReferenceData;
    private final PurpleConfig mConfig;

    private final String inputVCF;
    private boolean mEnabled;
    private final String mOutputVCF;
    private final TumorMutationalLoad mTumorMutationalLoad;
    private final MicrosatelliteIndels mMicrosatelliteIndels;
    private final SomaticVariantDrivers mDrivers;
    private final SomaticVariantFactory mSomaticVariantFactory;
    private final RChartData mRChartData;
    private final DriverGenePanel mGenePanel;
    private final List<PeakModel> mPeakModel;
    private final Set<String> mReportedGenes;
    private final List<VariantContext> mDownsampledSnp;
    private final List<VariantContext> mDownsampledIndel;
    private final int mSnpMod;
    private final int mIndelMod;
    private int mSnpCount;
    private int mIndelCount;

    public SomaticStream(
            final PurpleConfig config, final ReferenceData referenceData,
            int snpCount, int indelCount, final List<PeakModel> peakModel, final String somaticVcfFilename)
    {
        mReferenceData = referenceData;
        mConfig = config;

        mSnpMod = snpCount <= MAX_DOWNSAMPLE ? 1 : snpCount / MAX_DOWNSAMPLE;
        mIndelMod = indelCount <= MAX_DOWNSAMPLE ? 1 : indelCount / MAX_DOWNSAMPLE;

        mGenePanel = referenceData.GenePanel;
        mPeakModel = peakModel;
        mOutputVCF = config.OutputDir + config.TumorId + ".purple.somatic.vcf.gz";
        mEnabled = !somaticVcfFilename.isEmpty();
        inputVCF = somaticVcfFilename;
        mTumorMutationalLoad = new TumorMutationalLoad();
        mMicrosatelliteIndels = new MicrosatelliteIndels();
        mDrivers = new SomaticVariantDrivers(mGenePanel);
        mSomaticVariantFactory = SomaticVariantFactory.passOnlyInstance();
        mRChartData = new RChartData(config.OutputDir, config.TumorId);

        mReportedGenes = Sets.newHashSet();
        mDownsampledSnp = Lists.newArrayList();
        mDownsampledIndel = Lists.newArrayList();
    }

    public double microsatelliteIndelsPerMb()
    {
        return mMicrosatelliteIndels.microsatelliteIndelsPerMb();
    }

    @NotNull
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

    @NotNull
    public TumorMutationalStatus tumorMutationalBurdenPerMbStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromBurdenPerMb(tumorMutationalBurdenPerMb()) : TumorMutationalStatus.UNKNOWN;
    }

    @NotNull
    public TumorMutationalStatus tumorMutationalLoadStatus()
    {
        return mEnabled ? TumorMutationalStatus.fromLoad(tumorMutationalLoad()) : TumorMutationalStatus.UNKNOWN;
    }

    @NotNull
    public List<DriverCatalog> drivers(@NotNull final List<GeneCopyNumber> geneCopyNumbers)
    {
        return mDrivers.build(geneCopyNumbers);
    }

    @NotNull
    public Set<String> reportedGenes()
    {
        return mReportedGenes;
    }

    @NotNull
    public List<VariantContext> downsampledVariants()
    {
        List<VariantContext> result = Lists.newArrayList();
        result.addAll(mDownsampledIndel);
        result.addAll(mDownsampledSnp);

        return result;
    }

    public void processAndWrite(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions) throws IOException
    {
        final Consumer<VariantContext> downsampleConsumer = context ->
        {
            VariantContextDecorator variant = new VariantContextDecorator(context);
            if(variant.type() == VariantType.INDEL)
            {
                mIndelCount++;
                if(mIndelCount % mIndelMod == 0)
                {
                    mDownsampledIndel.add(context);
                }
            }
            else
            {
                mSnpCount++;
                if(mSnpCount % mSnpMod == 0)
                {
                    mDownsampledSnp.add(context);
                }
            }
        };

        final Consumer<VariantContext> driverConsumer =
                x -> mSomaticVariantFactory.createVariant(mConfig.TumorId, x).ifPresent(somatic ->
                {
                    boolean reported = mDrivers.add(somatic);
                    if(reported)
                    {
                        x.getCommonInfo().putAttribute(REPORTED_FLAG, true);
                        mReportedGenes.add(new VariantContextDecorator(x).gene());
                    }
                });

        if(mEnabled)
        {
            try (
                    VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false);
                    VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(mOutputVCF)
                            .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                            .build())
            {
                final Consumer<VariantContext> consumer = mTumorMutationalLoad.andThen(mMicrosatelliteIndels)
                        .andThen(driverConsumer)
                        .andThen(writer::add)
                        .andThen(mRChartData)
                        .andThen(downsampleConsumer);

                final SomaticVariantEnrichment enricher = new SomaticVariantEnrichment(
                        mConfig.DriverEnabled,
                        mConfig.SomaticFitting.clonalityBinWidth(),
                        mConfig.Version,
                        mConfig.ReferenceId,
                        mConfig.TumorId,
                        mReferenceData.RefGenome,
                        purityAdjuster,
                        mGenePanel,
                        copyNumbers,
                        fittedRegions,
                        mReferenceData.SomaticHotspots,
                        mReferenceData.TranscriptRegions,
                        mPeakModel,
                        consumer);

                final VCFHeader header = enricher.enrichHeader(vcfReader.getFileHeader());
                writer.writeHeader(header);

                for(VariantContext context : vcfReader)
                {
                    enricher.accept(context);
                }

                enricher.flush();
                mRChartData.write();
            }
        }

    }
}
