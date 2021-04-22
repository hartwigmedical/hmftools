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
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Downsample;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

class CircosCharts {

    private static final int MAX_PLOT_POINTS = 25000;

    private final ExecutorService executorService;
    private final String referenceSample;
    private final String tumorSample;
    private final ChartConfig config;
    private final String baseCircosTumorSample;
    private final String baseCircosReferenceSample;
    private final boolean isHg38;

    CircosCharts(final ConfigSupplier configSupplier, final ExecutorService executorService) {
        this.tumorSample = configSupplier.commonConfig().tumorSample();
        this.referenceSample = configSupplier.commonConfig().refSample();
        this.config = configSupplier.chartConfig();
        this.executorService = executorService;
        this.isHg38 = configSupplier.refGenomeConfig().isHg38();
        this.baseCircosTumorSample = config.circosDirectory() + File.separator + tumorSample;
        this.baseCircosReferenceSample = config.circosDirectory() + File.separator + referenceSample;
    }

    void write(@NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumber,
            @NotNull final List<VariantContext> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<FittedRegion> regions, @NotNull final List<AmberBAF> bafs) throws IOException {

        final List<VariantContextDecorator> somatics = somaticVariants.stream()
                .map(VariantContextDecorator::new)
                .filter(VariantContextDecorator::isPass)
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
    public List<Future<Integer>> chartFutures() {
        final List<Future<Integer>> futures = Lists.newArrayList();
        final Optional<String> circosBinary = config.circosBinary();
        if (circosBinary.isPresent()) {
            futures.add(executorService.submit(() -> generateCircos(circosBinary.get(), "input")));
            futures.add(executorService.submit(() -> generateCircos(circosBinary.get(), "circos")));
        }

        return futures;
    }

    @Nullable
    private Integer generateCircos(@NotNull final String executable, @NotNull final String type) throws IOException, InterruptedException {
        CircosExecution execution = new CircosExecution(executable);

        final String inputConfig = confFile(type);
        final String outputPath = config.plotDirectory();
        final String outputFile = tumorSample + "." + type + ".png";

        return execution.generateCircos(inputConfig, outputPath, outputFile, config.circosDirectory());
    }

    @NotNull
    private String confFile(@NotNull final String type) {
        return baseCircosTumorSample + "." + type + "." + "conf";
    }

    private void writeStructuralVariants(@NotNull final List<StructuralVariant> structuralVariants) throws IOException {
        CircosLinkWriter.writeVariants(baseCircosTumorSample + ".link.circos", structuralVariants);
    }

    private void writeFittedRegions(@NotNull final List<FittedRegion> fittedRegions) throws IOException {
        CircosFileWriter.writeRegions(baseCircosReferenceSample + ".ratio.circos",
                fittedRegions.stream().filter(x -> Doubles.positive(x.unnormalisedObservedNormalRatio())).collect(Collectors.toList()),
                ObservedRegion::unnormalisedObservedNormalRatio);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".ratio.circos",
                fittedRegions.stream().filter(x -> Doubles.positive(x.observedTumorRatio())).collect(Collectors.toList()),
                ObservedRegion::observedTumorRatio);
    }

    private void writeBafs(@NotNull final List<AmberBAF> bafs) throws IOException {
        CircosFileWriter.writePositions(baseCircosTumorSample + ".baf.circos", bafs, AmberBAF::tumorBAF);
    }

    private void writeCopyNumbers(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".map.circos", copyNumbers, x -> x.minorAlleleCopyNumber() - 1);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".cnv.circos", copyNumbers, x -> x.averageTumorCopyNumber() - 2);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".baf.circos", copyNumbers, PurpleCopyNumber::averageActualBAF);
    }

    private void writeEnrichedSomatics(@NotNull final List<VariantContextDecorator> somaticVariants) throws IOException {
        CircosSNPWriter.writePositions(baseCircosTumorSample + ".snp.circos", Downsample.downsample(MAX_PLOT_POINTS, snp(somaticVariants)));
        CircosINDELWriter.writePositions(baseCircosTumorSample + ".indel.circos",
                Downsample.downsample(MAX_PLOT_POINTS, indel(somaticVariants)));
    }

    private void writeConfig(@NotNull final Gender gender) throws IOException {
        writeConfig(gender, "circos");
        writeConfig(gender, "input");
        if (isHg38) {
            copyResourceToCircos("gaps.38.txt", "gaps.txt");
        } else {
            copyResourceToCircos("gaps.37.txt", "gaps.txt");
        }
    }

    private void writeConfig(@NotNull final Gender gender, @NotNull final String type) throws IOException {
        final Charset charset = StandardCharsets.UTF_8;
        String content = readResource("/circos/" + type + ".template");
        content = content.replaceAll("SAMPLE", tumorSample);
        content = content.replaceAll("REFERENCE", referenceSample);
        content = content.replaceAll("EXCLUDE", gender.equals(Gender.FEMALE) ? "hsY" : "hsZ");
        content = content.replaceAll("KARYOTYPE",
                isHg38 ? "data/karyotype/karyotype.human.hg38.txt" : "data/karyotype/karyotype.human.hg19.txt");
        Files.write(new File(confFile(type)).toPath(), content.getBytes(charset));
    }

    private void copyResourceToCircos(@NotNull final String inputName, @NotNull final String outputName) throws IOException {
        Charset charset = StandardCharsets.UTF_8;
        final String content = readResource("/circos/" + inputName);
        final String outputFilename = config.circosDirectory() + File.separator + outputName;
        Files.write(new File(outputFilename).toPath(), content.getBytes(charset));
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    @NotNull
    private List<VariantContextDecorator> snp(@NotNull final List<VariantContextDecorator> somaticVariants) {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.SNP)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }

    @NotNull
    private List<VariantContextDecorator> indel(@NotNull final List<VariantContextDecorator> somaticVariants) {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.INDEL)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }
}
