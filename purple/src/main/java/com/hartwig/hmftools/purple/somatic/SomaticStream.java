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
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;
import com.hartwig.hmftools.purple.plot.RChartData;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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

    private final SomaticFitConfig somaticFitConfig;
    private final CommonConfig commonConfig;
    private final String inputVCF;
    private boolean mEnabled;
    private final String outputVCF;
    private final TumorMutationalLoad tumorMutationalLoad;
    private final MicrosatelliteIndels microsatelliteIndels;
    private final SomaticVariantDrivers drivers;
    private final SomaticVariantFactory somaticVariantFactory;
    private final RChartData rChartData;
    private final DriverGenePanel genePanel;
    private final List<PeakModel> peakModel;
    private final Set<String> reportedGenes = Sets.newHashSet();
    private final List<VariantContext> downsampledSnp = Lists.newArrayList();
    private final List<VariantContext> downsampledIndel = Lists.newArrayList();
    private final int snpMod;
    private final int indelMod;
    private int snpCount;
    private int indelCount;

    public SomaticStream(
            final PurpleConfig config, final ReferenceData referenceData,
            int snpCount, int indelCount, final List<PeakModel> peakModel, final String somaticVcfFilename)
    {
        mReferenceData = referenceData;
        mConfig = config;

        this.snpMod = snpCount <= MAX_DOWNSAMPLE ? 1 : snpCount / MAX_DOWNSAMPLE;
        this.indelMod = indelCount <= MAX_DOWNSAMPLE ? 1 : indelCount / MAX_DOWNSAMPLE;

        this.genePanel = referenceData.GenePanel;
        this.somaticFitConfig = config.somaticConfig();
        this.commonConfig = config.commonConfig();
        this.peakModel = peakModel;
        this.outputVCF = commonConfig.outputDirectory() + File.separator + commonConfig.tumorSample() + ".purple.somatic.vcf.gz";
        this.mEnabled = !somaticVcfFilename.isEmpty();
        this.inputVCF = somaticVcfFilename;
        this.tumorMutationalLoad = new TumorMutationalLoad();
        this.microsatelliteIndels = new MicrosatelliteIndels();
        this.drivers = new SomaticVariantDrivers(genePanel);
        this.somaticVariantFactory = SomaticVariantFactory.passOnlyInstance();
        this.rChartData = new RChartData(commonConfig.outputDirectory(), commonConfig.tumorSample());
    }

    public double microsatelliteIndelsPerMb()
    {
        return microsatelliteIndels.microsatelliteIndelsPerMb();
    }

    @NotNull
    public MicrosatelliteStatus microsatelliteStatus()
    {
        return mEnabled ? MicrosatelliteStatus.fromIndelsPerMb(microsatelliteIndelsPerMb()) : MicrosatelliteStatus.UNKNOWN;
    }

    public double tumorMutationalBurdenPerMb()
    {
        return tumorMutationalLoad.burdenPerMb();
    }

    public int tumorMutationalLoad()
    {
        return tumorMutationalLoad.load();
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
        return drivers.build(geneCopyNumbers);
    }

    @NotNull
    public Set<String> reportedGenes()
    {
        return reportedGenes;
    }

    @NotNull
    public List<VariantContext> downsampledVariants()
    {
        List<VariantContext> result = Lists.newArrayList();
        result.addAll(downsampledIndel);
        result.addAll(downsampledSnp);

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
                indelCount++;
                if(indelCount % indelMod == 0)
                {
                    downsampledIndel.add(context);
                }
            }
            else
            {
                snpCount++;
                if(snpCount % snpMod == 0)
                {
                    downsampledSnp.add(context);
                }
            }
        };

        final Consumer<VariantContext> driverConsumer =
                x -> somaticVariantFactory.createVariant(commonConfig.tumorSample(), x).ifPresent(somatic ->
                {
                    boolean reported = drivers.add(somatic);
                    if(reported)
                    {
                        x.getCommonInfo().putAttribute(REPORTED_FLAG, true);
                        reportedGenes.add(new VariantContextDecorator(x).gene());
                    }
                });

        if(mEnabled)
        {
            try (
                    VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false);
                    VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                            .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                            .build())
            {

                final Consumer<VariantContext> consumer = tumorMutationalLoad.andThen(microsatelliteIndels)
                        .andThen(driverConsumer)
                        .andThen(writer::add)
                        .andThen(rChartData)
                        .andThen(downsampleConsumer);

                final SomaticVariantEnrichment enricher = new SomaticVariantEnrichment(
                        mConfig.DriverEnabled,
                        somaticFitConfig.clonalityBinWidth(),
                        commonConfig.version(),
                        commonConfig.tumorSample(),
                        mReferenceData.RefGenome,
                        purityAdjuster,
                        genePanel,
                        copyNumbers,
                        fittedRegions,
                        mReferenceData.SomaticHotspots,
                        mReferenceData.TranscriptRegions,
                        peakModel,
                        consumer);

                final VCFHeader header = enricher.enrichHeader(vcfReader.getFileHeader());
                writer.writeHeader(header);

                for(VariantContext context : vcfReader)
                {
                    enricher.accept(context);
                }

                enricher.flush();
                rChartData.write();
            }
        }

    }
}
