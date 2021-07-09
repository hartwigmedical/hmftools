package com.hartwig.hmftools.purple.plot;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.circos.CircosFileWriter;
import com.hartwig.hmftools.common.circos.CircosINDELWriter;
import com.hartwig.hmftools.common.circos.CircosLinkWriter;
import com.hartwig.hmftools.common.circos.CircosSNPWriter;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Downsample;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class CircosCharts
{
    private static final int MAX_PLOT_POINTS = 25000;

    private final ExecutorService mExecutorService;
    private final ChartConfig mConfig;
    private final boolean mIsHg38;

    private String mCurrentSampleId;
    private String mCurrentReferenceId;
    private String mBaseCircosTumorSample;
    private String mBaseCircosReferenceSample;

    public CircosCharts(final PurpleConfig configSupplier, final ExecutorService executorService, boolean isHg38)
    {
        mConfig = configSupplier.Charting;
        mExecutorService = executorService;
        mIsHg38 = isHg38;
    }

    public void write(
            final String referenceId, final String sampleId,
            final Gender gender, final List<PurpleCopyNumber> copyNumber,
            final List<VariantContextDecorator> somaticVariants, final List<StructuralVariant> structuralVariants,
            final List<FittedRegion> regions, final List<AmberBAF> bafs) throws IOException
    {
        mCurrentReferenceId = referenceId;
        mCurrentSampleId = sampleId;
        mBaseCircosTumorSample = mConfig.CircosDirectory + File.separator + mCurrentSampleId;
        mBaseCircosReferenceSample = mConfig.CircosDirectory + File.separator + mCurrentReferenceId;

        final List<VariantContextDecorator> somatics = somaticVariants.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .collect(Collectors.toList());

        writeConfig(gender);
        writeCopyNumbers(copyNumber);
        writeEnrichedSomatics(somatics);
        writeStructuralVariants(structuralVariants);
        writeFittedRegions(Downsample.downsample(MAX_PLOT_POINTS, regions));
        writeBafs(Downsample.downsample(MAX_PLOT_POINTS, bafs));
    }

    @NotNull
    public List<Future<Integer>> chartFutures()
    {
        final List<Future<Integer>> futures = Lists.newArrayList();
        final Optional<String> circosBinary = mConfig.CircosBinary;
        if(circosBinary.isPresent())
        {
            futures.add(mExecutorService.submit(() -> generateCircos(circosBinary.get(), "input")));
            futures.add(mExecutorService.submit(() -> generateCircos(circosBinary.get(), "circos")));
        }

        return futures;
    }

    @Nullable
    private Integer generateCircos(final String executable, final String type) throws IOException, InterruptedException
    {
        CircosExecution execution = new CircosExecution(executable);

        final String inputConfig = confFile(type);
        final String outputPath = mConfig.PlotDirectory;
        final String outputFile = mCurrentSampleId + "." + type + ".png";

        return execution.generateCircos(inputConfig, outputPath, outputFile, mConfig.CircosDirectory);
    }

    @NotNull
    private String confFile(final String type)
    {
        return mBaseCircosTumorSample + "." + type + "." + "conf";
    }

    private void writeStructuralVariants(final List<StructuralVariant> structuralVariants) throws IOException
    {
        CircosLinkWriter.writeVariants(mBaseCircosTumorSample + ".link.circos", structuralVariants);
    }

    private void writeFittedRegions(final List<FittedRegion> fittedRegions) throws IOException
    {
        CircosFileWriter.writeRegions(mBaseCircosReferenceSample + ".ratio.circos",
                fittedRegions.stream().filter(x -> Doubles.positive(x.unnormalisedObservedNormalRatio())).collect(Collectors.toList()),
                ObservedRegion::unnormalisedObservedNormalRatio);
        CircosFileWriter.writeRegions(mBaseCircosTumorSample + ".ratio.circos",
                fittedRegions.stream().filter(x -> Doubles.positive(x.observedTumorRatio())).collect(Collectors.toList()),
                ObservedRegion::observedTumorRatio);
    }

    private void writeBafs(final List<AmberBAF> bafs) throws IOException
    {
        CircosFileWriter.writePositions(mBaseCircosTumorSample + ".baf.circos", bafs, AmberBAF::tumorBAF);
    }

    private void writeCopyNumbers(final List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        CircosFileWriter.writeRegions(mBaseCircosTumorSample + ".map.circos", copyNumbers, x -> x.minorAlleleCopyNumber() - 1);
        CircosFileWriter.writeRegions(mBaseCircosTumorSample + ".cnv.circos", copyNumbers, x -> x.averageTumorCopyNumber() - 2);
        CircosFileWriter.writeRegions(mBaseCircosTumorSample + ".baf.circos", copyNumbers, PurpleCopyNumber::averageActualBAF);
    }

    private void writeEnrichedSomatics(final List<VariantContextDecorator> somaticVariants) throws IOException
    {
        CircosSNPWriter.writePositions(mBaseCircosTumorSample + ".snp.circos", Downsample.downsample(MAX_PLOT_POINTS, snp(somaticVariants)));
        CircosINDELWriter.writePositions(mBaseCircosTumorSample + ".indel.circos",
                Downsample.downsample(MAX_PLOT_POINTS, indel(somaticVariants)));
    }

    private void writeConfig(final Gender gender) throws IOException
    {
        writeConfig(gender, "circos");
        writeConfig(gender, "input");
        if(mIsHg38)
        {
            copyResourceToCircos("gaps.38.txt", "gaps.txt");
        }
        else
        {
            copyResourceToCircos("gaps.37.txt", "gaps.txt");
        }
    }

    private void writeConfig(final Gender gender, final String type) throws IOException
    {
        final Charset charset = StandardCharsets.UTF_8;
        String content = readResource("/circos/" + type + ".template");
        content = content.replaceAll("SAMPLE", mCurrentSampleId);
        content = content.replaceAll("REFERENCE", mCurrentReferenceId);
        content = content.replaceAll("EXCLUDE", gender.equals(Gender.FEMALE) ? "hsY" : "hsZ");
        content = content.replaceAll("KARYOTYPE",
                mIsHg38 ? "data/karyotype/karyotype.human.hg38.txt" : "data/karyotype/karyotype.human.hg19.txt");
        Files.write(new File(confFile(type)).toPath(), content.getBytes(charset));
    }

    private void copyResourceToCircos(final String inputName, final String outputName) throws IOException
    {
        Charset charset = StandardCharsets.UTF_8;
        final String content = readResource("/circos/" + inputName);
        final String outputFilename = mConfig.CircosDirectory + File.separator + outputName;
        Files.write(new File(outputFilename).toPath(), content.getBytes(charset));
    }

    @NotNull
    private String readResource(final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    @NotNull
    private List<VariantContextDecorator> snp(final List<VariantContextDecorator> somaticVariants)
    {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.SNP)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }

    @NotNull
    private List<VariantContextDecorator> indel(final List<VariantContextDecorator> somaticVariants)
    {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.INDEL)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }
}
