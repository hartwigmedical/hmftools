package com.hartwig.hmftools.patientreporter.viralbreakend;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViralBreakendFactory {

    private static final Logger LOGGER = LogManager.getLogger(ViralBreakendFactory.class);
    private static final String DELIMITER = "\t";

    private ViralBreakendFactory() {
    }

    @NotNull
    public static List<Viralbreakend> readViralBreakend(@NotNull String viralBreakendTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(viralBreakendTsv).toPath()));
    }

    @NotNull
    private static List<Viralbreakend> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(ViralBreakendFactory::fromString)
                .collect(toList());
    }

    @NotNull
    private static Viralbreakend fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutableViralbreakend.Builder builder = ImmutableViralbreakend.builder()
                .taxidGenus(values[0])
                .nameGenus(values[1])
                .readsGenusTree(values[2])
                .taxidSpecies(values[3])
                .nameSpecies(values[4])
                .readsSpeciesTree(values[5])
                .taxidAssigned(values[6])
                .nameAssigned(values[7])
                .readsAssignedTree(values[8])
                .readsAssignedDirect(values[9])
                .Reference(values[10])
                .referenceTaxid(values[11])
                .referenceKmerCount(values[12])
                .alternateKmerCountRname(values[13])
                .startpos(values[14])
                .endpos(values[15])
                .numreads(values[16])
                .covbases(values[17])
                .coverage(values[18])
                .meandepth(values[19])
                .meanbaseq(values[20])
                .meanmapq(values[21])
                .integrations(values[22])
                .QCStatus(values[23]);

        return builder.build();
    }
}