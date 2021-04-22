package com.hartwig.hmftools.cobalt.diploid;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

class DiploidCount implements Comparable<DiploidCount> {

    private static final String DELIMITER = "\t";

    @NotNull
    public static Map<GenomePosition, DiploidCount> readDiploidCountAsMap(final String inputFile) throws IOException {
        return Files.readAllLines(new File(inputFile).toPath())
                .stream()
                .map(DiploidCount::new)
                .collect(Collectors.toMap(x -> x.position, x -> x));
    }

    @NotNull
    public static List<DiploidCount> readDiploidCountAsList(final String inputFile) throws IOException {
        return Files.readAllLines(new File(inputFile).toPath()).stream().map(DiploidCount::new).collect(Collectors.toList());
    }

    private final GenomePosition position;
    private int diploid;
    private int count;

    private DiploidCount(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        this.position = GenomePositions.create(values[0], Long.parseLong(values[1]));
        this.diploid = Integer.parseInt(values[2]);
        this.count = Integer.parseInt(values[3]);
    }

    DiploidCount(final GenomePosition position, final int diploid, final int count) {
        this.position = position;
        this.diploid = diploid;
        this.count = count;
    }

    void incrementTotal() {
        count++;
    }

    void incrementDiploid() {
        diploid++;
    }

    public int getDiploid() {
        return diploid;
    }

    public int getCount() {
        return count;
    }

    public String chromosome() {
        return position.chromosome();
    }

    public long position() {
        return position.position();
    }

    public double proportionIsDiploid(int count) {
        return 1d * getDiploid() / count;
    }

    @NotNull
    @Override
    public String toString() {
        return new StringJoiner(DELIMITER).add(position.chromosome())
                .add(String.valueOf(position.position()))
                .add(String.valueOf(diploid))
                .add(String.valueOf(count))
                .toString();
    }

    @Override
    public int compareTo(@NotNull final DiploidCount o) {
        return position.compareTo(o.position);
    }
}
