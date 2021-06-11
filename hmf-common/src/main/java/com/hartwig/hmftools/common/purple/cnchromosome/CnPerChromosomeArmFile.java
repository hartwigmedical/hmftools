package com.hartwig.hmftools.common.purple.cnchromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class CnPerChromosomeArmFile {

    private CnPerChromosomeArmFile() {

    }

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION = ".cnv.chromosomearm.somatic.tsv";
    private static final String DELIMITER = "\t";

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + COPYNUMBER_PER_CHROMOSOME_ARM_EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull List<CnPerChromosomeArmData> cnPerChromosomeArm) throws IOException {
        Files.write(new File(filename).toPath(), toLines(cnPerChromosomeArm));
    }

    @VisibleForTesting
    @NotNull
    public static List<String> toLines(@NotNull final List<CnPerChromosomeArmData> cnPerChromosomeArm) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());

        for (CnPerChromosomeArmData chromsomeArm : cnPerChromosomeArm) {
            lines.add(new StringJoiner(DELIMITER).add(chromsomeArm.chromosome().toString())
                    .add(chromsomeArm.chromosomeArm().name())
                    .add(FORMAT.format(chromsomeArm.copyNumber()))
                    .toString());
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome").add("chromosomeArm").add("copyNumber").toString();
    }
}
