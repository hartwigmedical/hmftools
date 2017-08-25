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
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.circos.CircosFileWriter;
import com.hartwig.hmftools.common.circos.CircosLinkWriter;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.CircosConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.CollectionUtil;

class GenerateCircosDataHelper {

    private static final Logger LOGGER = LogManager.getLogger(GenerateCircosDataHelper.class);

    private static final int MAX_SOMATIC_VARIANTS = 40000;
    private final String sample;
    private final CircosConfig config;
    private final String baseCircosSample;
    private final String circosConfigOutput;

    GenerateCircosDataHelper(final String sample, final CircosConfig config) throws IOException {
        this.sample = sample;
        this.config = config;
        this.baseCircosSample = config.circosDirectory() + File.separator + sample;
        this.circosConfigOutput = baseCircosSample + ".circos.conf";
    }

    private void createDirectory(final String dir) throws IOException {
        final File output = new File(dir);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create circos output directory " + dir);
        }
    }

    void write(@NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumber,
            @NotNull final List<EnrichedSomaticVariant> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants)
            throws IOException, InterruptedException {

        createDirectory(config.plotDirectory());
        createDirectory(config.circosDirectory());

        writeConfig(gender);
        writeCopyNumbers(copyNumber);
        writeEnrichedSomatics(somaticVariants);
        writeStructualVariants(structuralVariants);

        if (config.circosBinary().isPresent()) {
            generateCircos(config.circosBinary().get());
        }
    }

    private void generateCircos(@NotNull final String executable) throws IOException, InterruptedException {
        final File errorFile = new File(baseCircosSample + ".circos.error");
        final File outputFile = new File(baseCircosSample + ".circos.out");
        final String plotFileName = sample + ".circos.png";

        final String[] command = new String[8];
        command[0] = executable;
        command[1] = "-nosvg";
        command[2] = "-conf";
        command[3] = new File(circosConfigOutput).getAbsolutePath();
        command[4] = "-outputdir";
        command[5] = new File(config.plotDirectory()).getAbsolutePath();
        command[6] = "-outputfile";
        command[7] = plotFileName;

        LOGGER.info(String.format("Generating CIRCOS via command: %s", CollectionUtil.join(Arrays.asList(command), " ")));
        int result = new ProcessBuilder(command).redirectError(errorFile).redirectOutput(outputFile).start().waitFor();
        if (result != 0) {
            LOGGER.warn("Error creating circos plot. Examine error file {} for details.", errorFile.toString());
        }
    }

    private void writeStructualVariants(@NotNull final List<StructuralVariant> structuralVariants) throws IOException {
        CircosLinkWriter.writeVariants(baseCircosSample + ".link.circos", structuralVariants);
    }

    private void writeCopyNumbers(@NotNull final List<PurpleCopyNumber> copyNumber) throws IOException {
        CircosFileWriter.writeRegions(baseCircosSample + ".cnv.circos", copyNumber, x -> x.averageTumorCopyNumber() - 2);
        CircosFileWriter.writeRegions(baseCircosSample + ".baf.circos", copyNumber, PurpleCopyNumber::averageActualBAF);
    }

    private void writeEnrichedSomatics(@NotNull final List<EnrichedSomaticVariant> somaticVariants) throws IOException {
        final List<EnrichedSomaticVariant> downsampledSomaticVariants = downsample(filter(somaticVariants));
        CircosFileWriter.writePositions(baseCircosSample + ".snp.circos", downsampledSomaticVariants, EnrichedSomaticVariant::adjustedVAF);
    }

    private void writeConfig(@NotNull final Gender gender) throws IOException {
        Charset charset = StandardCharsets.UTF_8;
        String content = readResource("/circos/circos.template");
        content = content.replaceAll("SAMPLE", sample);
        content = content.replaceAll("EXCLUDE", gender.equals(Gender.MALE) ? "hsZ" : "hsY");
        Files.write(new File(circosConfigOutput).toPath(), content.getBytes(charset));

        copyResourceToCircos("gaps.txt");
        copyResourceToCircos("ideogram.conf");
        copyResourceToCircos("ticks.conf");
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
    private List<EnrichedSomaticVariant> filter(@NotNull final List<EnrichedSomaticVariant> somaticVariants) {
        return somaticVariants.stream()
                .filter(x -> x.type() == VariantType.SNP)
                .filter(x -> HumanChromosome.fromString(x.chromosome()).intValue() <= 25)
                .collect(Collectors.toList());
    }

    private List<EnrichedSomaticVariant> downsample(List<EnrichedSomaticVariant> variants) {
        if (variants.size() <= MAX_SOMATIC_VARIANTS) {
            return variants;
        }

        long scale = Math.round(Math.ceil(1.0 * variants.size() / MAX_SOMATIC_VARIANTS));
        final List<EnrichedSomaticVariant> result = Lists.newArrayList();
        for (int i = 0; i < variants.size(); i++) {
            if (i % scale == 0) {
                result.add(variants.get(i));
            }
        }

        return result;
    }
}
