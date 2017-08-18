package com.hartwig.hmftools.purple;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.circos.CircosFileWriter;
import com.hartwig.hmftools.common.circos.CircosLinkWriter;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

class GenerateCircosDataHelper {

    private static final int MAX_SOMATIC_VARIANTS = 40000;
    private final String baseDataOutput;
    private final String outputDirectory;
    private final String sample;

    GenerateCircosDataHelper(final String sample, final String outputDirectory) throws IOException {
        this.sample = sample;
        this.outputDirectory = outputDirectory;
        this.baseDataOutput = outputDirectory + File.separator + sample;
        final File output = new File(outputDirectory);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create circos output directory " + outputDirectory);
        }
    }

    void write(@NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumber, @NotNull final List<EnrichedSomaticVariant> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants)
            throws IOException {
        writeConfig(gender);
        writeCopyNumbers(copyNumber);
        writeEnrichedSomatics(somaticVariants);
        writeStructualVariants(structuralVariants);
    }

    private void writeStructualVariants(@NotNull final List<StructuralVariant> structuralVariants) throws IOException {
        CircosLinkWriter.writeVariants(baseDataOutput + ".link.circos", structuralVariants);
    }

    private void writeCopyNumbers(@NotNull final List<PurpleCopyNumber> copyNumber) throws IOException {
        CircosFileWriter.writeRegions(baseDataOutput + ".cnv.circos", copyNumber, x -> x.averageTumorCopyNumber() - 2);
        CircosFileWriter.writeRegions(baseDataOutput + ".baf.circos", copyNumber, PurpleCopyNumber::averageActualBAF);
    }

    private void writeEnrichedSomatics(@NotNull final List<EnrichedSomaticVariant> somaticVariants) throws IOException {
        final List<EnrichedSomaticVariant> downsampledSomaticVariants = downsample(filter(somaticVariants));
        CircosFileWriter.writePositions(baseDataOutput + ".snp.circos", downsampledSomaticVariants, EnrichedSomaticVariant::adjustedVAF);
    }

    private void writeConfig(@NotNull final Gender gender) throws IOException {
        final String circosConfigOutput = baseDataOutput + ".circos.conf";
        Charset charset = StandardCharsets.UTF_8;
        String content = readResource("/circos/circos.template");
        content = content.replaceAll("SAMPLE", sample);
        content = content.replaceAll("EXCLUDE", gender.equals(Gender.MALE) ? "hsZ" : "hsY");
        Files.write(new File(circosConfigOutput).toPath(), content.getBytes(charset));

        copyResourceToFile("gaps.txt");
        copyResourceToFile("ideogram.conf");
        copyResourceToFile("ticks.conf");
    }

    private void copyResourceToFile(@NotNull final String name) throws IOException {
        Charset charset = StandardCharsets.UTF_8;
        final String content = readResource("/circos/" + name);
        final String outputFilename = outputDirectory + File.separator + name;
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
                .filter(x -> Chromosomes.asInt(x.chromosome()) <= 25)
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
