package com.hartwig.hmftools.virusinterpreter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class TaxonomyDbFile {

    private static final Logger LOGGER = LogManager.getLogger(TaxonomyDbFile.class);

    private static final String SEPARATOR = "\t";

    private TaxonomyDbFile() {
    }

    @NotNull
    public static TaxonomyDb loadFromTsv(@NotNull String taxonomyDbTsv) throws IOException {
        List<String> linesTaxonomyDb = Files.readAllLines(new File(taxonomyDbTsv).toPath());

        Map<Integer, String> taxidToNameMap = Maps.newHashMap();

        for (String line : linesTaxonomyDb) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                int taxid = Integer.parseInt(parts[0].trim());
                String name = parts[1].trim();
                taxidToNameMap.put(taxid, name);
            } else {
                LOGGER.warn("Suspicious line detected in taxonomy db tsv: {}", line);
            }
        }

        return new TaxonomyDb(taxidToNameMap);
    }
}
