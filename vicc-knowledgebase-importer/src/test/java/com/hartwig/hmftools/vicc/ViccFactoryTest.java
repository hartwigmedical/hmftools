package com.hartwig.hmftools.vicc;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class ViccFactoryTest {

    private static final Logger LOGGER = LogManager.getLogger(ViccFactoryTest.class);

    @Test
    @Ignore
    public void canReadJsonKnowledgeBaseFile() throws IOException {
        final String baseDir =
                System.getProperty("user.home") + File.separator + "hmf" + File.separator + "projects" + File.separator + "vicc";
        final String inputFile = baseDir + File.separator + "all.json";

        List<ViccEntry> viccEntries = ViccFactory.readViccKnowledgebaseJsonFile(inputFile);

        assertNotNull(viccEntries);

        analyzeViccEntries(viccEntries);
    }

    private static void analyzeViccEntries(@NotNull List<ViccEntry> viccEntries) {
        LOGGER.info("Analyzing VICC entries - total count = " + viccEntries.size());

        LOGGER.info("Analyzing entries per SOURCE");
        countPerCategory(viccEntries, entry -> Lists.newArrayList(entry.source()));

        LOGGER.info("Analyzing entries per TAG");
        countPerCategory(viccEntries, ViccEntry::tags);

        LOGGER.info("Analyzing entries per DEVTAG");
        countPerCategory(viccEntries, ViccEntry::devTags);

        LOGGER.info("Analyzing entries per GENE");
        countPerCategory(viccEntries, ViccEntry::genes);
    }

    private static void countPerCategory(@NotNull List<ViccEntry> viccEntries,
            @NotNull Function<ViccEntry, List<String>> keyGenerator) {
        Map<String, Integer> countsPerKey = Maps.newHashMap();
        for (ViccEntry entry : viccEntries) {
            List<String> keys = keyGenerator.apply(entry);
            for (String key : keys) {
                countsPerKey.merge(key, 1, (a, b) -> a + b);
            }
        }

        Set<Integer> sortedCounts = Sets.newTreeSet(Comparator.reverseOrder());
        sortedCounts.addAll(countsPerKey.values());

        for (Integer count : sortedCounts) {
            for (Map.Entry<String, Integer> entry : countsPerKey.entrySet()) {
                if (entry.getValue().equals(count)) {
                    LOGGER.info("  " + entry.getKey() + ": " + entry.getValue() + " entries");
                }
            }
        }
    }
}