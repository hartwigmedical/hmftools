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
    public static List<EventType> determineEventType(@NotNull ViccEntry viccEntry) {

        boolean combinedEvent = false;
        String biomarkerType = Strings.EMPTY;
        String gene = Strings.EMPTY;
        String eventSource = Strings.EMPTY;
        List<String> event = Lists.newArrayList();

        List<EventType> eventType = Lists.newArrayList();

        Source type = Source.sourceFromKnowledgebase(viccEntry.source());

        for (Feature feature : viccEntry.features()) {
            switch (type) {
                case ONCOKB: // extract info oncokb
                    gene = feature.geneSymbol();
                    eventSource = feature.name();
                    biomarkerType = feature.biomarkerType();
                    if (eventSource.contains("+")) {
                        String[] combinedEventConvertToSingleEvent = eventSource.split(" ", 2);
                        if (combinedEventConvertToSingleEvent.length == 2) {
                            combinedEvent = true;
                            event.add(combinedEventConvertToSingleEvent[0]);
                            event.add(combinedEventConvertToSingleEvent[1]);
                        }
                    } else {
                        event.add(eventSource);
                    }

                    //    LOGGER.info(event);

                    //                    if (event.equals("NA")) {
                    //                        String[] eventArray = feature.name().split(" ", 2);
                    //                        if (eventArray.length == 1) {
                    //                            event = eventArray[0];
                    //                        } else if (eventArray.length == 2) {
                    //                            if (eventArray[1].equals("Fusion")) {
                    //                                event = eventArray[1];
                    //                            } else {
                    //                                event = feature.name();
                    //                            }
                    //                        }
                    //                        if (event.contains("Exon")) {
                    //                            event = feature.name();
                    //                        } else if (Pattern.compile("[0-9]").matcher(event).find()) {
                    //                            event = "manual curated mutation";
                    //                        }
                    //                    }
                    //
                    //                    if (feature.provenanceRule() == null) {
                    //                        LOGGER.info("No provencence rule known");
                    //                    } else if (feature.provenanceRule().equals("gene_only")) {
                    //                        event = "gene_level";
                    //                    }
                    //
                    //                    if (event.isEmpty()) {
                    //                        LOGGER.warn(
                    //                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                                feature.name(),
                    //                                feature.geneSymbol(),
                    //                                feature.biomarkerType(),
                    //                                feature.description(),
                    //                                type,
                    //                                event);
                    //                    }

                    break;
                case CGI: // extract info for cgi
                    gene = feature.geneSymbol();
                    eventSource = feature.name();
                    biomarkerType = feature.biomarkerType();
                    if (eventSource.contains("+")) {
                        String[] combinedEventConvertToSingleEvent = eventSource.split(" ", 2);
                        if (combinedEventConvertToSingleEvent.length == 2) {
                            combinedEvent = true;
                            event.add(combinedEventConvertToSingleEvent[0]);
                            event.add(combinedEventConvertToSingleEvent[1]);
                        }
                    } else {
                        event.add(eventSource);
                    }
                    LOGGER.info(event);

                    //
                    //                    if (feature.provenanceRule() == null) {
                    //                        LOGGER.info("No provencence rule known");
                    //                    } else if (feature.provenanceRule().equals("gene_only")) {
                    //                        event = "gene_level";
                    //                    }
                    //
                    //                    if (event.isEmpty()) {
                    //                        LOGGER.warn(
                    //                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                                feature.name(),
                    //                                feature.geneSymbol(),
                    //                                feature.biomarkerType(),
                    //                                feature.description(),
                    //                                type,
                    //                                event);
                    //                    }
                    break;
                case CIVIC: // extract info for civic
                    gene = feature.geneSymbol();
                    eventSource = feature.name();
                    biomarkerType = feature.biomarkerType();
                    if (eventSource.contains("+")) {
                        String[] combinedEventConvertToSingleEvent = eventSource.split(" ", 2);
                        if (combinedEventConvertToSingleEvent.length == 2) {
                            combinedEvent = true;
                            event.add(combinedEventConvertToSingleEvent[0]);
                            event.add(combinedEventConvertToSingleEvent[1]);
                        }
                    } else {
                        event.add(eventSource);
                    }
                    //                    if (event == null || event.equals("N/A")) {
                    //                        event = feature.name();
                    //                    }
                    //
                    //                    if (feature.provenanceRule() == null) {
                    //                        LOGGER.info("No provencence rule known");
                    //                    } else if (feature.provenanceRule().equals("gene_only")) {
                    //                        event = "gene_level";
                    //                    }
                    //
                    //                    if (event.contains("Prime") || event.contains("Exon") || event.contains("EXON") || event.contains("Frameshift") || event
                    //                            .contains("EXPRESSION") || event.equals("RS34743033") || event.contains("SPLICE VARIANT") || event.contains(
                    //                            "PHOSPHORYLATION")) {
                    //                        event = event;
                    //                    } else if (Pattern.compile("[0-9]").matcher(event).find()) {
                    //                        event = "manual curated mutation";
                    //                    }
                    //
                    //                    if (event.isEmpty()) {
                    //                        LOGGER.warn(
                    //                                "Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' and description {} on {} and event is {}",
                    //                                feature.name(),
                    //                                feature.geneSymbol(),
                    //                                feature.biomarkerType(),
                    //                                feature.description(),
                    //                                type,
                    //                                event);
                    //                    }

                    break;
                case BRCA: // extract info for brca
                    break;
                case JAX: // extract info for jax
                    break;
                case JAX_TRIALS: // extract info for jax trials
                    break;
                case MOLECULARMATCH: // extract info for molecular match
                    break;
                case MOLECULARMATCH_TRIALS: // extract info for molecular match trials
                    break;
                case PMKB: // extract info for pmkb
                    break;
                case SAGE: // extract info for sage
                    break;
                default:
                    LOGGER.warn("Unknown knowledgebase!");
            }

            eventType.add(ImmutableEventType.builder()
                    .gene(gene)
                    .combinedEvent(combinedEvent)
                    .event(eventSource)
                    .interpretEventType(event)
                    .biomarkerType(biomarkerType).build());
        }
        return eventType;
    }
}
