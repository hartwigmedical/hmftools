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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Downsample;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.ChartConfig;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.Nullable;

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
            final List<ObservedRegion> regions, final List<AmberBAF> bafs) throws IOException
    {
        mCurrentReferenceId = referenceId;
        mCurrentSampleId = sampleId;
        mBaseCircosTumorSample = mConfig.CircosDirectory + mCurrentSampleId;

        mBaseCircosReferenceSample = mCurrentReferenceId != null ? mConfig.CircosDirectory + mCurrentReferenceId : "";

        final List<VariantContextDecorator> somatics = somaticVariants.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .collect(Collectors.toList());

        writeConfig(gender);
        writeCopyNumbers(copyNumber);
        writeEnrichedSomatics(somatics);
        writeStructuralVariants(structuralVariants);
        writeObservedRegions(Downsample.downsample(MAX_PLOT_POINTS, regions));
        writeBafs(Downsample.downsample(MAX_PLOT_POINTS, bafs));
    }

    public List<Future<Integer>> chartFutures()
    {
        final List<Future<Integer>> futures = Lists.newArrayList();
        futures.add(mExecutorService.submit(() -> generateCircos(mConfig.CircosBinary, "input")));
        futures.add(mExecutorService.submit(() -> generateCircos(mConfig.CircosBinary, "circos")));

        return futures;
    }

    @Nullable
    private Integer generateCircos(final String executable, final String type) throws IOException, InterruptedException
    {
        CircosExecution execution = new CircosExecution(executable);

        final String inputConfig = confFile(type);
        final String outputPath = mConfig.PlotDirectory;
        final String outputFile = mCurrentSampleId + "." + type + ".png";

        return execution.generateCircos(inputConfig, outputPath, outputFile);
    }

    private String confFile(final String type)
    {
        return mBaseCircosTumorSample + "." + type + "." + "conf";
    }

    private void writeStructuralVariants(final List<StructuralVariant> structuralVariants) throws IOException
    {
        CircosLinkWriter.writeVariants(mBaseCircosTumorSample + ".link.circos", structuralVariants);
    }

    private void writeObservedRegions(final List<ObservedRegion> fittedRegions) throws IOException
    {
        List<ObservedRegion> selectRegions = fittedRegions.stream()
                .filter(x -> x.germlineStatus() == GermlineStatus.DIPLOID).collect(Collectors.toList());

        // make the glyph size proportional to number of cobalt windows
        ToDoubleFunction<ObservedRegion> glyphSizeFunc = (ObservedRegion region) -> {
            int windowCount = region.depthWindowCount();
            return 4 * (Math.log10(Math.max(windowCount, 1)) + 1);
        };

        if(!mBaseCircosReferenceSample.isEmpty())
        {
            CircosFileWriter.writeRegions(mBaseCircosReferenceSample + ".ratio.circos",
                    selectRegions.stream().filter(x -> Doubles.positive(x.unnormalisedObservedNormalRatio())).collect(Collectors.toList()),
                    ObservedRegion::unnormalisedObservedNormalRatio,
                    glyphSizeFunc);
        }

        CircosFileWriter.writeRegions(
                mBaseCircosTumorSample + ".ratio.circos",
                selectRegions.stream().filter(x -> Doubles.positive(x.observedTumorRatio())).collect(Collectors.toList()),
                ObservedRegion::observedTumorRatio,
                glyphSizeFunc);
    }

    private void writeBafs(final List<AmberBAF> bafs) throws IOException
    {
        CircosFileWriter.writePositions(mBaseCircosTumorSample + ".baf.circos", bafs, x -> x.TumorBAF);
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
        content = content.replaceAll("REFERENCE", mCurrentReferenceId != null ? mCurrentReferenceId : mCurrentSampleId);

        // use the tumor colour for reference in tumor-only mode to keep those points the same
        content = content.replaceAll("COBALT_REF_COLOUR", mCurrentReferenceId != null ? "green" : "cobalt");

        content = content.replaceAll("EXCLUDE", gender.equals(Gender.FEMALE) ? "hsY" : "hsZ");
        content = content.replaceAll("KARYOTYPE",
                mIsHg38 ? "data/karyotype/karyotype.human.hg38.txt" : "data/karyotype/karyotype.human.hg19.txt");
        Files.write(new File(confFile(type)).toPath(), content.getBytes(charset));
    }

    private void copyResourceToCircos(final String inputName, final String outputName) throws IOException
    {
        Charset charset = StandardCharsets.UTF_8;
        final String content = readResource("/circos/" + inputName);
        final String outputFilename = mConfig.CircosDirectory + outputName;
        Files.write(new File(outputFilename).toPath(), content.getBytes(charset));
    }

    private String readResource(final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    private List<VariantContextDecorator> snp(final List<VariantContextDecorator> somaticVariants)
    {
        return somaticVariants.stream().filter(x -> x.type() == VariantType.SNP).collect(Collectors.toList());
    }

    private List<VariantContextDecorator> indel(final List<VariantContextDecorator> somaticVariants)
    {
        return somaticVariants.stream().filter(x -> x.type() == VariantType.INDEL).collect(Collectors.toList());
    }
}
