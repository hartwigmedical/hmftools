package com.hartwig.hmftools.knowledgebasegenerator.cnv;

import java.util.List;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class GeneratingCNV {

    private static final Logger LOGGER = LogManager.getLogger(GeneratingCNV.class);

    private static final String DELIMITER = "\t";
    private static final String SOURCE_LINK_SEPARATOR = ";";
    private static final String NEW_LINE = "\n";

    private static final List<String> AMPLIFICATION =
            Lists.newArrayList("amplification", "Amplification", "Gain", "overexpression", "amp", "over exp");
    private static final List<String> DELETION =
            Lists.newArrayList("deletion", "Deletion", "Copy Number Loss", "Loss", "loss", "undexpression");

    public static void generatingCNVs(@NotNull ViccEntry viccEntries) {
        String headerActionableCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type"
                        + DELIMITER + "Cancer Type" + DELIMITER + "Level" + DELIMITER + "Direction" + NEW_LINE;
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;
        extractCNV(viccEntries);
    }

    public static void extractCNV(@NotNull ViccEntry viccEntry) {
        String event = Strings.EMPTY;
        for (Feature feature : viccEntry.features()) {
            if (viccEntry.source().equals("oncokb")) { // extract info oncokb
                event = feature.biomarkerType();
                if (event.equals("NA")) {
                    event = feature.description().substring(feature.description().lastIndexOf(" ") + 1);
                    if (Pattern.compile("[0-9]").matcher(event).find()) {
                        LOGGER.info("hand curated: " + event);
                        event = "mutation";
                    }
                }
                LOGGER.info("event oncokb: " + event);
            } else if (viccEntry.source().equals("cgi")) { // extract info for cgi
                event = feature.biomarkerType();
                LOGGER.info("event cgi: " + event);


            }

        }

    }

}
