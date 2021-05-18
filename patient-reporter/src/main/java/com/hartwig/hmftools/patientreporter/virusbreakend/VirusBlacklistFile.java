package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusBlacklistFile {

    private static final Logger LOGGER = LogManager.getLogger(TaxonomyDbFile.class);

    private static final String SEPARATOR = "\t";

    private VirusBlacklistFile() {
    }

    @NotNull
    public static VirusBlacklistModel buildFromTsv(@NotNull String virusBlacklistTsv) throws IOException {
        List<String> linesVirusBlacklist = Files.readAllLines(new File(virusBlacklistTsv).toPath());

        Set<Integer> blacklistedGenera = Sets.newHashSet();
        Set<Integer> blacklistedSpecies = Sets.newHashSet();

        for (String line : linesVirusBlacklist.subList(1, linesVirusBlacklist.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 3) {
                int taxid = Integer.parseInt(parts[0].trim());
                String taxidType = parts[1].trim();
                if (taxidType.equals("taxid_genus")) {
                    blacklistedGenera.add(taxid);
                } else if (taxidType.equals("taxid_species")) {
                    blacklistedSpecies.add(taxid);
                } else {
                    LOGGER.warn("Could not interpret virus blacklist entry of type '{}'", taxidType);
                }
            } else {
                LOGGER.warn("Suspicious line detected in virus blacklist tsv: {}", line);
            }
        }

        return new VirusBlacklistModel(blacklistedGenera, blacklistedSpecies);
    }
}
