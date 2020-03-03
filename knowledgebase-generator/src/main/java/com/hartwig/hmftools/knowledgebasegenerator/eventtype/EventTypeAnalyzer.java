package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class EventTypeAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);

    @NotNull
    public static List<EventType> determineEventType(@NotNull ViccEntry viccEntry, int num) {
        String event = Strings.EMPTY;
        List<EventType> eventType = Lists.newArrayList();
        String gene = Strings.EMPTY;
        String genomicMutation = Strings.EMPTY;

        Source type = Source.sourceFromKnowledgebase(viccEntry.source());

        for (Feature feature : viccEntry.features()) {
            switch (type) {
                case ONCOKB: // extract info oncokb
                    gene = feature.geneSymbol();
                    genomicMutation = feature.name();
                    event = feature.biomarkerType();

                    if (event.equals("NA")) {
                        String[] eventArray = feature.name().split(" ", 2);
                        if (eventArray.length == 1) {
                            event = eventArray[0];
                        } else if (eventArray.length == 2) {
                            if (eventArray[1].equals("Fusion")) {
                                event = eventArray[1];
                            } else {
                                event = feature.name();
                            }
                        }
                        if (event.contains("Exon")) {
                            event = feature.name();
                        } else if (Pattern.compile("[0-9]").matcher(event).find()) {
                            event = "manual curated mutation";
                        }
                    }

                    if (feature.provenanceRule() == null) {
                        LOGGER.info("No provencence rule known");
                    } else if (feature.provenanceRule().equals("gene_only")) {
                        event = "gene_level";
                    }

                    if (event.isEmpty()) {
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                                feature.name(),
                                feature.geneSymbol(),
                                feature.biomarkerType(),
                                feature.description(),
                                type,
                                event);
                    }

                    break;
                case CGI: // extract info for cgi
                    gene = feature.geneSymbol();
                    genomicMutation = feature.name();
                    event = feature.biomarkerType();

                    if (feature.provenanceRule() == null) {
                        LOGGER.info("No provencence rule known");
                    } else if (feature.provenanceRule().equals("gene_only")) {
                        event = "gene_level";
                    }

                    if (event.isEmpty()) {
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                                feature.name(),
                                feature.geneSymbol(),
                                feature.biomarkerType(),
                                feature.description(),
                                type,
                                event);
                    }
                    break;
                case BRCA: // extract info for brca //TODO
                    event = Strings.EMPTY;
                    //    LOGGER.info(num + ": event brca: " + event);
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    break;
                case CIVIC: // extract info for civic
                    gene = feature.geneSymbol();
                    genomicMutation = feature.name();
                    event = feature.biomarkerType();
                    if (event == null || event.equals("N/A")) {
                        event = feature.name();
                    }

                    if (feature.provenanceRule() == null) {
                        LOGGER.info("No provencence rule known");
                    } else if (feature.provenanceRule().equals("gene_only")) {
                        event = "gene_level";
                    }

                    if (event.contains("Prime") || event.contains("Exon") || event.contains("EXON") || event.contains("Frameshift") || event
                            .contains("EXPRESSION") || event.equals("RS34743033") || event.contains("SPLICE VARIANT") || event.contains(
                            "PHOSPHORYLATION")) {
                        event = event;
                    } else if (Pattern.compile("[0-9]").matcher(event).find()) {
                        event = "manual curated mutation";
                    }

                    if (event.isEmpty()) {
                        LOGGER.warn(
                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                                feature.name(),
                                feature.geneSymbol(),
                                feature.biomarkerType(),
                                feature.description(),
                                type,
                                event);
                    }

                    break;
                case JAX: // extract info for jax //TODO
                    event = feature.biomarkerType();
                    if (event == null) {
                        if (feature.description() != null) {
                            event = feature.description();
                            if (Pattern.compile("[0-9]").matcher(event).find()) {
                                event = event + " : description manual curated mutation";
                            }
                        } else {
                            event = feature.name();
                            if (Pattern.compile("[0-9]").matcher(event).find()) {
                                event = event + " : name manual curated mutation";
                            }
                        }
                    }
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    //  LOGGER.info(num + ": event jax: " + event);
                    break;
                case JAX_TRIALS: // extract info for jax trials //TODO
                    event = feature.biomarkerType();
                    if (event == null) {
                        if (feature.description() != null) {
                            event = feature.description();
                            if (Pattern.compile("[0-9]").matcher(event).find()) {
                                event = event + " :  manual curated mutation";
                            }
                        } else {
                            event = feature.name();
                            if (Pattern.compile("[0-9]").matcher(event).find()) {
                                event = event + " : name manual curated mutation";
                            }
                        }

                    }
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    //   LOGGER.info(num + ": event jax trials: " + event);
                    break;
                case MOLECULARMATCH: // extract info for molecular match //TODO

                    //                    LOGGER.info(num + "skip because it is from source molecularmatch");
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    break;
                case MOLECULARMATCH_TRIALS: // extract info for molecular match trials //TODO
                    event = feature.biomarkerType();
                    if (event == null) {
                        String[] eventArray = feature.description().split(" ");
                        if (eventArray.length == 1) {
                            event = event + "array: " + feature.description().split(" ", 2)[0];
                        } else {
                            event = event + "array: " + feature.description().split(" ", 2)[1];
                        }
                    }
                    //                    LOGGER.info(num + ": event molecular match trials: " + event);
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    break;
                case PMKB: // extract info for pmkb //TODO
                    event = feature.biomarkerType();
                    //   LOGGER.info(num + ": event pmkb: " + event);
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    break;
                case SAGE: // extract info for sage //TODO
                    event = feature.biomarkerType();
                    if (event == null) {
                        event = feature.description();
                    }
                    // LOGGER.info(num + ": event sage: " + event);
                    //                    LOGGER.info(
                    //                            "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                            feature.name(),
                    //                            feature.geneSymbol(),
                    //                            feature.biomarkerType(),
                    //                            feature.description(),
                    //                            type,
                    //                            event);
                    break;
                default:
                    LOGGER.warn(num + ": Unknown knowledgebase");
            }
            if (event == null) {
                event = Strings.EMPTY;
            }
            eventType.add(ImmutableEventType.builder().gene(gene).genomicMutation(genomicMutation).eventType(event).build());
        }
        return eventType;
    }
}
