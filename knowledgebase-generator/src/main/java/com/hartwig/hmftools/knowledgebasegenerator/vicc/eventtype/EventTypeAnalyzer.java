package com.hartwig.hmftools.knowledgebasegenerator.vicc.eventtype;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.tools.StringUtils;

public final class EventTypeAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);
    private static final String ONCOGENIC_MUTATION = "oncogenic mutation";
    private static final String FUSION_PAIR = "fusion pair";
    private static final String FUSION_PROMISCUOUS = "fusion promiscuous";

    private EventTypeAnalyzer() {
    }

    @NotNull
    public static List<EventType> determineEventTypes(@NotNull ViccEntry viccEntry) {
        boolean combinedEvent = false;
        String biomarkerType = Strings.EMPTY;
        String gene = Strings.EMPTY;
        String name = Strings.EMPTY;
        String description = Strings.EMPTY;
        String eventInfo;

        List<EventType> eventType = Lists.newArrayList();

        for (Feature feature : viccEntry.features()) {
            Map<String, List<String>> eventMap = Maps.newHashMap();

            switch (viccEntry.source()) {
                case ONCOKB: // extract info oncokb
                    name = feature.name();
                    gene = feature.geneSymbol();
                    biomarkerType = feature.biomarkerType();
                    description = feature.description();

                    if (name.equals("Fusions")) {
                        name = FUSION_PROMISCUOUS;
                    } else if (name.contains("Fusion")) {
                        if (name.contains(" - ")) {
                            gene = name.split(" Fusion")[0];
                            name = FUSION_PAIR;
                        } else {
                            gene = name.split(" ")[0];
                            name = FUSION_PAIR;
                        }
                    }

                    eventMap.put(gene, Lists.newArrayList(name));

                    if (eventMap.isEmpty()) {
                        LOGGER.warn("Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}' ",
                                feature.name(),
                                gene,
                                biomarkerType,
                                viccEntry.source());
                    }

                    break;
                case CGI: // extract info for cgi
                    name = feature.name();
                    biomarkerType = feature.biomarkerType();
                    description = feature.description();
                    if (name.contains("+")) {

                        String[] combinedEventConvertToSingleEvent = name.split(" \\+ ", 2);
                        gene = combinedEventConvertToSingleEvent[0].split(" ", 2)[0];
                        eventInfo = combinedEventConvertToSingleEvent[0].split(" ", 2)[1];
                        if (eventInfo.equals(".")) {
                            eventInfo = ONCOGENIC_MUTATION;
                        }

                        String geneCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[0];
                        String eventInfoCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[1];

                        //I assume, a combined event for actionability has 2 events. If more events, this will be not interpretated
                        if (combinedEventConvertToSingleEvent.length == 2) {
                            combinedEvent = true;

                            if (eventMap.size() == 0) {
                                eventMap.put(gene, Lists.newArrayList(eventInfo));
                                if (eventMap.containsKey(geneCombined)) {
                                    eventMap.put(geneCombined, Lists.newArrayList(eventInfo, eventInfoCombined));
                                } else {
                                    eventMap.put(gene, Lists.newArrayList(eventInfo));
                                    eventMap.put(geneCombined, Lists.newArrayList(eventInfoCombined));
                                }
                            }
                        } else if (combinedEventConvertToSingleEvent.length >= 2) {
                            LOGGER.warn("This event has more events, which is not interpretated!");
                        }
                    } else {
                        if (name.contains(":")) {
                            gene = name.split(":", 2)[0];
                            eventInfo = name.split(":", 2)[1];
                            if (eventInfo.equals(".")) {
                                eventInfo = ONCOGENIC_MUTATION;
                            }
                            eventMap.put(gene, Lists.newArrayList(eventInfo));
                        } else {
                            gene = name.split(" ", 2)[0];

                            if (name.split(" ", 2).length == 1 && gene.contains("-") && !biomarkerType.equals("mutant")) {
                                eventInfo = FUSION_PAIR;
                                eventMap.put(gene, Lists.newArrayList(eventInfo));
                            } else {
                                if (name.contains("fusion") || name.contains("Fusion")) {
                                    if (gene.contains("-")) {
                                        eventInfo = FUSION_PAIR;
                                    } else {
                                        eventInfo = FUSION_PROMISCUOUS;
                                    }
                                } else if (name.split(" ", 2).length == 2) {
                                    if (name.split(" ", 2)[1].equals(".")) {
                                        eventInfo = ONCOGENIC_MUTATION;
                                    } else {
                                        eventInfo = name.split(" ", 2)[1];
                                    }
                                } else {
                                    eventInfo = name;
                                }
                                eventMap.put(gene, Lists.newArrayList(eventInfo));

                            }
                        }
                    }
                    if (eventMap.isEmpty()) {
                        LOGGER.warn("Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType,
                                viccEntry.source());
                    }

                    break;
                case CIVIC: // extract info for civic
                    name = feature.name();
                    biomarkerType = feature.biomarkerType() == null ? "null" : feature.biomarkerType();
                    gene = feature.geneSymbol();
                    description = feature.description();

                    // assume there are no gene with "-" in their name
                    int count = StringUtils.countMatches(name, "-");
                    if (count >= 2 && !name.equals("LOSS-OF-FUNCTION") && !name.equals("Gain-of-Function")) {
                        LOGGER.warn("Fix for gene '{}'", name);
                    }
                    if (!name.contains("DEL") && !name.contains("Splicing alteration") && !name.contains("EXON") && !name.contains("c.")
                            && !name.contains("MUT") && !name.equals("LOSS-OF-FUNCTION") && !name.equals("Gain-of-Function")
                            && !name.contains("C.") && !name.equals("N-TERMINAL FRAME SHIFT") && !name.equals("COPY-NEUTRAL LOSS OF HETEROZYGOSITY")) {

                        if (name.contains("-") && name.contains(" ")) {
                            String[] combinedEventConvertToSingleEvent = name.split(" ", 2);

                            String fusion = combinedEventConvertToSingleEvent[0];
                            String variant = combinedEventConvertToSingleEvent[1];
                            String geneVariant = fusion.split("-")[1];

                            //I assume, a combined event for actionability has 2 events. If more events, this will be not interpretated
                            if (combinedEventConvertToSingleEvent.length == 2) {
                                combinedEvent = true;

                                if (eventMap.size() == 0) {
                                    eventMap.put(fusion, Lists.newArrayList(FUSION_PAIR));
                                    if (eventMap.containsKey(geneVariant)) {
                                        eventMap.put(geneVariant, Lists.newArrayList(FUSION_PAIR, variant));
                                    } else {
                                        eventMap.put(fusion, Lists.newArrayList(FUSION_PAIR));
                                        eventMap.put(geneVariant, Lists.newArrayList(variant));
                                    }
                                }
                            } else if (combinedEventConvertToSingleEvent.length >= 2) {
                                LOGGER.warn("This event has more events, which is not interpretated!");
                            }
                        } else if (name.contains("-") && !biomarkerType.equals("Missense Variant")) {
                            eventMap.put(name, Lists.newArrayList(FUSION_PAIR));
                        } else if (name.equals("TRUNCATING FUSION")) {
                            eventMap.put(gene, Lists.newArrayList(name));
                        } else if (name.contains("FUSION") || name.contains("FUSIONS")) {
                            eventMap.put(gene, Lists.newArrayList(FUSION_PROMISCUOUS));
                        } else {
                            if (name.contains("+")) {

                                combinedEvent = true;
                                String[] combinedEventConvertToSingleEvent = name.replace("+", " ").split(" ", 2);

                                String event1 = combinedEventConvertToSingleEvent[0];
                                String event2 = combinedEventConvertToSingleEvent[1];

                                if (eventMap.size() == 0) {
                                    eventMap.put(gene, Lists.newArrayList(event1));
                                    if (eventMap.containsKey(gene)) {
                                        eventMap.put(gene, Lists.newArrayList(event1, event2));
                                    }
                                } else {
                                    eventMap.put(gene, Lists.newArrayList(name));
                                }

                            } else {
                                eventMap.put(gene, Lists.newArrayList(name));
                            }
                        }
                    } else if (name.contains("+") && !name.contains("c.") && !name.contains("C.")) {
                        combinedEvent = true;
                        String[] combinedEventConvertToSingleEvent = name.split("\\+", 2);
                        String event1 = combinedEventConvertToSingleEvent[0];
                        String event2 = combinedEventConvertToSingleEvent[1];

                        eventMap.put(gene, Lists.newArrayList(event1, event2));

                    } else {
                        eventMap.put(gene, Lists.newArrayList(name));
                    }

                    if (eventMap.isEmpty()) {
                        LOGGER.warn("Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}' on source '{}'",
                                feature.name(),
                                gene,
                                biomarkerType,
                                viccEntry.source());
                    }

                    break;
                case BRCA: // extract info for brca
                    break;
                case JAX: // extract info for jax
                    break;
                case JAX_TRIALS: // extract info for jax trials
                    break;
                case MOLECULAR_MATCH: // extract info for molecular match
                    break;
                case MOLECULAR_MATCH_TRIALS: // extract info for molecular match trials
                    break;
                case PMKB: // extract info for pmkb
                    break;
                case SAGE: // extract info for sage
                    break;
                default:
                    LOGGER.warn("Unknown knowledgebase!");
            }

            eventType.add(ImmutableEventType.builder()
                    .combinedEvent(combinedEvent)
                    .biomarkerType(biomarkerType)
                    .gene(gene)
                    .name(name)
                    .description(description)
                    .source(viccEntry.source())
                    .eventMap(eventMap)
                    .build());
        }
        return eventType;
    }
}
