package com.hartwig.hmftools.common.virus;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class VirusBreakendFile {

    private static final String DELIMITER = "\t";

    private VirusBreakendFile() {
    }

    @NotNull
    public static List<VirusBreakend> read(@NotNull String virusBreakendTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(virusBreakendTsv).toPath()));
    }

    @NotNull
    private static List<VirusBreakend> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(VirusBreakendFile::fromString)
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
                .RName(values[14])
                .startPos(Integer.parseInt(values[15]))
                .endPos(Integer.parseInt(values[16]))
                .numReads(Integer.parseInt(values[17]))
                .covBases(Integer.parseInt(values[18]))
                .coverage(Double.parseDouble(values[19]))
                .meanDepth(Double.parseDouble(values[20]))
                .meanBaseQ(Double.parseDouble(values[21]))
                .meanMapQ(Double.parseDouble(values[22]))
                .integrations(Integer.parseInt(values[23]))
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES);

        if (values.length == 25) {
            builder.qcStatus(VirusBreakendQCStatus.convert(values[24]));
        }

        return builder.build();
    }
}