package com.hartwig.hmftools.protect.variants.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GermlineReportingFile {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingFile.class);
    private static final String SEPARATOR = "\t";

    private GermlineReportingFile() {
    }

    @NotNull
    public static GermlineReportingModel buildFromTsv(@NotNull String germlineGenesCsv) throws IOException {
        List<String> linesGermlineGenes = Files.readAllLines(new File(germlineGenesCsv).toPath());

            Map<String, GermlineReporting> germlineGenesAndNotifyMap = Maps.newHashMap();

            for (String line : linesGermlineGenes.subList(1, linesGermlineGenes.size())) {
                String[] parts = line.split(SEPARATOR);

                if (parts.length == 3) {
                    String gene = parts[0].trim();
                    String notifyGene = parts[1].trim().toLowerCase();
                    ConditionReportingVariant conditon = ConditionReportingVariant.fromConditionString(parts[2].trim());
                    String variant = Strings.EMPTY;
                    germlineGenesAndNotifyMap =
                            extractGermlineInfo(germlineGenesAndNotifyMap, notifyGene, gene, conditon, variant);

                } else if (parts.length == 4) {
                    String gene = parts[0].trim();
                    String notifyGene = parts[1].trim().toLowerCase();
                    ConditionReportingVariant conditon = ConditionReportingVariant.fromConditionString(parts[2].trim());
                    String variant = parts[3].trim();
                    germlineGenesAndNotifyMap =
                            extractGermlineInfo(germlineGenesAndNotifyMap, notifyGene, gene, conditon, variant);
                } else {
                    LOGGER.warn("Suspicious line detected in germline reporting genes: {}", line);
                }
            }

            return new GermlineReportingModel(germlineGenesAndNotifyMap);
        }

        @NotNull
        public static Map<String, GermlineReporting> extractGermlineInfo(@NotNull Map<String, GermlineReporting> germlineGenesAndNotifyMap,
                @NotNull String notifyGene, @NotNull String gene, @NotNull ConditionReportingVariant condition,
                @NotNull String variant) {
            GermlineReporting reporting;

            switch (notifyGene) {
                case "true":
                    reporting = ImmutableGermlineReporting.builder()
                            .notifyClinicalGeneticus(true)
                            .condition(condition)
                            .variant(variant)
                            .build();

                    germlineGenesAndNotifyMap.put(gene, reporting);
                    break;
                case "false":
                    reporting = ImmutableGermlineReporting.builder()
                            .notifyClinicalGeneticus(false)
                            .condition(condition)
                            .variant(variant)
                            .build();
                    germlineGenesAndNotifyMap.put(gene, reporting);
                    break;
                default:
                    LOGGER.warn("Could not interpret notification string in germline reporting genes: {}", notifyGene);
            }
            return germlineGenesAndNotifyMap;
        }
    }


