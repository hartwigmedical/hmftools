package com.hartwig.hmftools.purple;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.baf.TumorBAF;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.circos.CircosFileWriter;
import com.hartwig.hmftools.common.circos.CircosLinkWriter;
import com.hartwig.hmftools.common.circos.CircosSNPWriter;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.CircosConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.CollectionUtil;

class GenerateCircosData {

    private static final Logger LOGGER = LogManager.getLogger(GenerateCircosData.class);

    private static final int MAX_PLOT_POINTS = 25000;
    private final String referenceSample;
    private final String tumorSample;
    private final CircosConfig config;
    private final String baseCircosTumorSample;
    private final String baseCircosReferenceSample;

    GenerateCircosData(final ConfigSupplier configSupplier) throws IOException {
        this.tumorSample = configSupplier.commonConfig().tumorSample();
        this.referenceSample = configSupplier.commonConfig().refSample();
        this.config = configSupplier.circosConfig();
        this.baseCircosTumorSample = config.circosDirectory() + File.separator + tumorSample;
        this.baseCircosReferenceSample = config.circosDirectory() + File.separator + referenceSample;
    }

    private void createDirectory(final String dir) throws IOException {
        final File output = new File(dir);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create circos output directory " + dir);
        }
    }

    void write(@NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumber,
            @NotNull final List<PurityAdjustedSomaticVariant> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<FittedRegion> regions, @NotNull final List<TumorBAF> bafs) throws IOException, InterruptedException {

        createDirectory(config.plotDirectory());
        createDirectory(config.circosDirectory());

        writeConfig(gender);
        writeCopyNumbers(copyNumber);
        writeEnrichedSomatics(somaticVariants);
        writeStructualVariants(structuralVariants);
        writeFittedRegions(downsample(regions));
        writeBafs(downsample(bafs));

        if (config.circosBinary().isPresent()) {
            generateCircos(config.circosBinary().get(), "ratio");
            generateCircos(config.circosBinary().get(), "circos");
        }
    }

    private void generateCircos(@NotNull final String executable, @NotNull final String type) throws IOException, InterruptedException {
        final File errorFile = file(type, "error");
        final File outputFile = file(type, "out");
        final String plotFileName = tumorSample + "." + type + ".png";

        final String[] command = new String[8];
        command[0] = executable;
        command[1] = "-nosvg";
        command[2] = "-conf";
        command[3] = confFile(type).getAbsolutePath();
        command[4] = "-outputdir";
        command[5] = new File(config.plotDirectory()).getAbsolutePath();
        command[6] = "-outputfile";
        command[7] = plotFileName;

        LOGGER.info(String.format("Generating " + type.toUpperCase() + " via command: %s",
                CollectionUtil.join(Arrays.asList(command), " ")));
        int result = new ProcessBuilder(command).redirectError(errorFile).redirectOutput(outputFile).start().waitFor();
        if (result != 0) {
            LOGGER.warn("Error creating circos plot. Examine error file {} for details.", errorFile.toString());
        }
    }

    @NotNull
    private File confFile(@NotNull final String type) {
        return file(type, "conf");
    }

    @NotNull
    private File file(@NotNull final String type, @NotNull final String extension) {
        return new File(baseCircosTumorSample + "." + type + "." + extension);
    }

    private void writeStructualVariants(@NotNull final List<StructuralVariant> structuralVariants) throws IOException {
        CircosLinkWriter.writeVariants(baseCircosTumorSample + ".link.circos", structuralVariants);
    }

    private void writeFittedRegions(@NotNull final List<FittedRegion> fittedRegions) throws IOException {
        CircosFileWriter.writeRegions(baseCircosReferenceSample + ".ratio.circos", fittedRegions, ObservedRegion::observedNormalRatio);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".ratio.circos", fittedRegions, ObservedRegion::observedTumorRatio);
    }
    private void writeBafs(@NotNull final List<TumorBAF> bafs) throws IOException {
        CircosFileWriter.writePositions(baseCircosTumorSample + ".baf.circos", bafs, TumorBAF::baf);
    }

    private void writeCopyNumbers(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".map.circos", copyNumbers, x -> minorAllelePloidy(x) - 1);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".cnv.circos", copyNumbers, x -> x.averageTumorCopyNumber() - 2);
        CircosFileWriter.writeRegions(baseCircosTumorSample + ".baf.circos", copyNumbers, PurpleCopyNumber::averageActualBAF);
    }

    private void writeEnrichedSomatics(@NotNull final List<PurityAdjustedSomaticVariant> somaticVariants) throws IOException {
        final List<PurityAdjustedSomaticVariant> downsampledSomaticVariants = downsample(filter(somaticVariants));
        CircosSNPWriter.writePositions(baseCircosTumorSample + ".snp.circos", downsampledSomaticVariants);
    }

    private double minorAllelePloidy(@NotNull final PurpleCopyNumber copyNumber) {
        return Doubles.lessThan(copyNumber.averageActualBAF(), 0.50)
                ? 0
                : Math.max(0, (1 - copyNumber.averageActualBAF()) * copyNumber.averageTumorCopyNumber());
    }

    private void writeConfig(@NotNull final Gender gender) throws IOException {
        writeConfig(gender, "circos");
        writeConfig(gender, "ratio");

        copyResourceToCircos("gaps.txt");
        copyResourceToCircos("ideogram.conf");
        copyResourceToCircos("ticks.conf");
    }

    private void writeConfig(@NotNull final Gender gender, @NotNull final String type) throws IOException {
        final Charset charset = StandardCharsets.UTF_8;
        String content = readResource("/circos/" + type + ".template");
        content = content.replaceAll("SAMPLE", tumorSample);
        content = content.replaceAll("REFERENCE", referenceSample);
        content = content.replaceAll("EXCLUDE", gender.equals(Gender.MALE) ? "hsZ" : "hsY");
        Files.write(confFile(type).toPath(), content.getBytes(charset));
    }

    private void copyResourceToCircos(@NotNull final String name) throws IOException {
        Charset charset = StandardCharsets.UTF_8;
        final String content = readResource("/circos/" + name);
        final String outputFilename = config.circosDirectory() + File.separator + name;
        Files.write(new File(outputFilename).toPath(), content.getBytes(charset));
    }

    private String readResource(@NotNull final String resource) throws IOException {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    @NotNull
    private List<PurityAdjustedSomaticVariant> filter(@NotNull final List<PurityAdjustedSomaticVariant> somaticVariants) {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.SNP)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }

    private static <T> List<T> downsample(List<T> variants) {
        if (variants.size() <= MAX_PLOT_POINTS) {
            return variants;
        }

        long scale = Math.round(Math.ceil(1.0 * variants.size() / MAX_PLOT_POINTS));
        final List<T> result = Lists.newArrayList();
        for (int i = 0; i < variants.size(); i++) {
            if (i % scale == 0) {
                result.add(variants.get(i));
            }
        }
        return result;
    }
}
