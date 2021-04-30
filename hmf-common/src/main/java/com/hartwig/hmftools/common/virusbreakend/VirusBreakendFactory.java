package com.hartwig.hmftools.common.virusbreakend;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class VirusBreakendFactory {

    private static final String DELIMITER = "\t";

    private VirusBreakendFactory() {
    }

    @NotNull
    public static List<VirusBreakend> readVirusBreakend(@NotNull String virusBreakendTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(virusBreakendTsv).toPath()));
    }

    @NotNull
    private static List<VirusBreakend> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(VirusBreakendFactory::fromString)
                .collect(toList());
    }

    @NotNull
    private static VirusBreakend fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutableVirusBreakend.Builder builder = ImmutableVirusBreakend.builder()
                .taxidGenus(Integer.parseInt(values[0]))
                .nameGenus(values[1])
                .readsGenusTree(Integer.parseInt(values[2]))
                .taxidSpecies(Integer.parseInt(values[3]))
                .nameSpecies(values[4])
                .readsSpeciesTree(Integer.parseInt(values[5]))
                .taxidAssigned(Integer.parseInt(values[6]))
                .nameAssigned(values[7])
                .readsAssignedTree(Integer.parseInt(values[8]))
                .readsAssignedDirect(Integer.parseInt(values[9]))
                .reference(values[10])
                .referenceTaxid(Integer.parseInt(values[11]))
                .referenceKmerCount(Integer.parseInt(values[12]))
                .alternateKmerCountRname(values[13])
                .startpos(Integer.parseInt(values[14]))
                .endpos(Integer.parseInt(values[15]))
                .numreads(Integer.parseInt(values[16]))
                .covbases(Integer.parseInt(values[17]))
                .coverage(Integer.parseInt(values[18]))
                .meandepth(Integer.parseInt(values[19]))
                .meanbaseq(Integer.parseInt(values[20]))
                .meanmapq(Integer.parseInt(values[21]))
                .integrations(Integer.parseInt(values[22]))
                .QCStatus(Integer.parseInt(values[23]));

        return builder.build();
    }
}