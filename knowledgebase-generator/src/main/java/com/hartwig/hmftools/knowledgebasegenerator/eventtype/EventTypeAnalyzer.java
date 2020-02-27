package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
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
        String sourceKnowledgebase = viccEntry.source();
        for (Feature feature : viccEntry.features()) {
            switch (sourceKnowledgebase) {
                case "oncokb": // extract info oncokb
                    OncoKb kbOncoKb = (OncoKb) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event.equals("NA")) {
                        event = feature.description().split(" ",2)[1];
                        if (Pattern.compile("[0-9]").matcher(event).find()) {
                          //  LOGGER.info("manual curated: " + event + " to mutation");
                            event = "manual curated mutation";
                        }
                    }
                    LOGGER.info(num + ": event oncokb: " + event);
                    break;
                case "cgi": // extract info for cgi
                    Cgi kbCgi = (Cgi) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event cgi: " + event);
                    break;
                case "brca": // extract info for brca  //TODO
                    Brca kbBrca = (Brca) kbSpecificObject;
                    event = Strings.EMPTY;
                      LOGGER.info(num + ": event brca: " + event);
                    break;
                case "civic": // extract info for civic //TODO
                    Civic kbCivic = (Civic) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        event = "manual curated mutation";
                    } else if (event.equals("N/A")) {
                        event = feature.description().split(" ",2)[1];
                    }

                    LOGGER.info(num + ": event civic: " + event);
                    break;
                case "jax": // extract info for jax //TODO
                    Jax kbJax = (Jax) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event jax: " + event);
                    break;
                case "jax_trials": // extract info for jax trials //TODO
                    JaxTrials kbJaxTrials = (JaxTrials) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event jax trials: " + event);
                    break;
                case "molecularmatch": // extract info for molecular match //TODO
                    MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;
                    event = feature.biomarkerType();
                    if (event == null) {
                        event = feature.description().split(" ",2)[1];
                    }
                    LOGGER.info(num + ": event molecular match: " + event);
                    break;
                case "molecularmatch_trials": // extract info for molecular match trials //TODO
                    MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event molecular match trials: " + event);
                    break;
                case "pmkb": // extract info for pmkb //TODO
                    Pmkb kbPmkb = (Pmkb) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event pmkb: " + event);
                    break;
                case "sage": // extract info for sage //TODO
                    Sage kbSage = (Sage) kbSpecificObject;
                    event = feature.biomarkerType();
                    LOGGER.info(num + ": event sage: " + event);
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
