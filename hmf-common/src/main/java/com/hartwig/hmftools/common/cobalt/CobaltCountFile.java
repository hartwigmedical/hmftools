package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.cobalt.CobaltRatioFile.fromLine;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public class CobaltCountFile {

    @NotNull
    public static Multimap<Chromosome, CobaltCount> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static Multimap<Chromosome, CobaltCount> fromLines(@NotNull final List<String> lines) {

        final Multimap<Chromosome, CobaltCount> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith("Ch")) {
                final CobaltRatio position = fromLine(line);
                result.put(HumanChromosome.fromString(position.chromosome()), position);
            }
        }

        return result;
    }

}
