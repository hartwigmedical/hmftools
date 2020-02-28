package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.knowledgebasegenerator.sourceKnowledgebase.Sources;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.brca.Brca;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class EventTypeAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);

    @NotNull
    public static List<EventType> determineEventType(@NotNull ViccEntry viccEntry, int num) {
        KbSpecificObject kbSpecificObject = viccEntry.KbSpecificObject();
        String event = Strings.EMPTY;
        List<EventType> eventType = Lists.newArrayList();

        Sources type = Sources.courceFromKnowledgebase(viccEntry.source());

        for (Feature feature : viccEntry.features()) {
            switch (type) {
                case ONCOKB: // extract info oncokb
                    OncoKb kbOncoKb = (OncoKb) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event.equals("NA")) {
                        String [] eventArray = feature.description().split(" ");
                        if (eventArray.length == 1) {
                            event = "array: " + feature.description().split(" ", 2)[0];
                        } else {
                            event = "array: " + feature.description().split(" ", 2)[1];
                        }

                        if (Pattern.compile("[0-9]").matcher(event).find()) {
                            event = "manual curated mutation";
                        }
                    }
                 //   LOGGER.info(num + ": event oncokb: " + event);
                    break;
                case CGI: // extract info for cgi
                    Cgi kbCgi = (Cgi) kbSpecificObject;
                    event = feature.biomarkerType();
                 //   LOGGER.info(num + ": event cgi: " + event);
                    break;
                case BRCA: // extract info for brca  //TODO
                    Brca kbBrca = (Brca) kbSpecificObject;
                    event = Strings.EMPTY;
                  //    LOGGER.info(num + ": event brca: " + event);
                    break;
                case CIVIC: // extract info for civic
                    Civic kbCivic = (Civic) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        event = "manual curated mutation";
                    } else if (event.equals("N/A")) {
                        String [] eventArray = feature.description().split(" ");
                        if (eventArray.length == 1) {
                            event = "array: " + feature.description().split(" ", 2)[0];
                        } else {
                            event = "array: " + feature.description().split(" ", 2)[1];
                            if (Pattern.compile("[0-9]").matcher(event).find()) {
                                event = feature.description().split(" ", 2)[1].split(" ", 2)[1];
                            }
                        }
                    }

                  //  LOGGER.info(num + ": event civic: " + event);
                    break;
                case JAX: // extract info for jax
                    Jax kbJax = (Jax) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        if (feature.description() != null){
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
                  //  LOGGER.info(num + ": event jax: " + event);
                    break;
                case JAX_TRIALS: // extract info for jax trials
                    JaxTrials kbJaxTrials = (JaxTrials) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        if (feature.description() != null){
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
                 //   LOGGER.info(num + ": event jax trials: " + event);
                    break;
                case MOLECULARMATCH: // extract info for molecular match //TODO what to do when features is empty
                    MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;

                    LOGGER.info(num + "skip because it is from source molecularmatch");
                    break;
                case MOLECULARMATCH_TRIALS: // extract info for molecular match trials //TODO
                    MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        String [] eventArray = feature.description().split(" ");
                        if (eventArray.length == 1) {
                            event = event + "array: " + feature.description().split(" ", 2)[0];
                        } else {
                            event = event + "array: " + feature.description().split(" ", 2)[1];
                        }
                    }
                    LOGGER.info(num + ": event molecular match trials: " + event);
                    break;
                case PMKB: // extract info for pmkb
                    Pmkb kbPmkb = (Pmkb) kbSpecificObject;
                    event = feature.biomarkerType();
                 //   LOGGER.info(num + ": event pmkb: " + event);
                    break;
                case SAGE: // extract info for sage
                    Sage kbSage = (Sage) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        event = feature.description();
                    }
                   // LOGGER.info(num + ": event sage: " + event);
                    break;
                default:
                    LOGGER.warn(num + ": Unknown knowledgebase");
            }
            if (event == null) {
                event = Strings.EMPTY;
            }
            eventType.add(ImmutableEventType.builder().gene("").eventType(event).build());
        }
        return eventType;
    }
}
