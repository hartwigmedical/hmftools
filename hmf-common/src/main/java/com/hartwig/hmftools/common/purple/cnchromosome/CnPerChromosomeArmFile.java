package com.hartwig.hmftools.common.purple.cnchromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class CnPerChromosomeArmFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION = ".cnv.chromosomearm.somatic.tsv";
    private static final String DELIMITER = "\t";

    private CnPerChromosomeArmFile() {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION;
    }

    public static void write(@NotNull String filename, @NotNull List<CnPerChromosomeArmData> cnPerChromosomeArm) throws IOException {
        Files.write(new File(filename).toPath(), toLines(cnPerChromosomeArm));
    }

    @NotNull
    private static List<String> toLines(@NotNull List<CnPerChromosomeArmData> cnPerChromosomeArm) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());

        for (CnPerChromosomeArmData chromosomeArm : cnPerChromosomeArm) {
            lines.add(new StringJoiner(DELIMITER).add(chromosomeArm.chromosome().toString())
                    .add(chromosomeArm.chromosomeArm().name())
                    .add(FORMAT.format(chromosomeArm.copyNumber()))
                    .toString());
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("chromosome").add("chromosomeArm").add("copyNumber").toString();
    }
}
