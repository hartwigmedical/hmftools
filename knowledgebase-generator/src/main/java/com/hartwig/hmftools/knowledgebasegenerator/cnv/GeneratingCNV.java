package com.hartwig.hmftools.knowledgebasegenerator.cnv;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
                "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type" + DELIMITER + "Cancer Type"
                        + DELIMITER + "Level" + DELIMITER + "Direction"  + NEW_LINE;
        String headerknownCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;
        extractCNV(viccEntries);
    }

    public static void extractCNV(@NotNull ViccEntry viccEntry) {

    }

}
