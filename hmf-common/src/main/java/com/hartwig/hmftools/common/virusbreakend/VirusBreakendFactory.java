package com.hartwig.hmftools.common.virusbreakend;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.logging.log4j.util.Strings;
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
        ImmutableVirusBreakend.Builder builder = ImmutableVirusBreakend.builder()
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
                .alternateKmerCount(Integer.parseInt(values[13]))
                .Rname(values[14])
                .startpos(Integer.parseInt(values[15]))
                .endpos(Integer.parseInt(values[16]))
                .numreads(Integer.parseInt(values[17]))
                .covbases(Integer.parseInt(values[18]))
                .coverage(Double.parseDouble(values[19]))
                .meandepth(Double.parseDouble(values[20]))
                .meanbaseq(Double.parseDouble(values[21]))
                .meanmapq(Double.parseDouble(values[22]))
                .integrations(Integer.parseInt(values[23]))
                .QCStatus(Strings.EMPTY);

        if (values.length == 25) {
            builder.QCStatus(values[24]);
        }

        return builder.build();
    }
}