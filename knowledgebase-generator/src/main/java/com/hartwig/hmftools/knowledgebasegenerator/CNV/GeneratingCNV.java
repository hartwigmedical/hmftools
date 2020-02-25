package com.hartwig.hmftools.knowledgebasegenerator.CNV;

import java.util.List;

import com.google.common.collect.Lists;
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
import org.jetbrains.annotations.NotNull;

public class GeneratingCNV {

    private static final Logger LOGGER = LogManager.getLogger(GeneratingCNV.class);

    private static final String DELIMTER = "\t";
    private static final String NEW_LINE = "\n";

    private static final List<String> AMPLIFICATION = Lists.newArrayList("amplification", "Amplification", "Gain", "overexpression", "amp", "over exp");
    private static final List<String> DELETION = Lists.newArrayList("deletion", "Deletion", "Copy Number Loss", "Loss", "loss", "undexpression");


    @NotNull
    public static void generatingCNVs(@NotNull List<ViccEntry> viccEntries) {

        String headerActionableCNV =
                "Gene" + DELIMTER + "Type" + DELIMTER + "Source" + DELIMTER + "Drug" + DELIMTER + "Drug Type" + DELIMTER + "Cancer Type"
                        + DELIMTER + "Level" + DELIMTER + "Direction" + DELIMTER + "Link" + NEW_LINE;
        String headerknownCNV =
                "Gene" + DELIMTER + "Type" + DELIMTER + "Source" + DELIMTER + "Drug" + DELIMTER + "Drug Type" + DELIMTER + "Cancer Type"
                        + DELIMTER + "Level" + DELIMTER + "Direction" + DELIMTER + "Link" + NEW_LINE;

        for (ViccEntry viccEntry : viccEntries) {
            KbSpecificObject kbSpecificObject = viccEntry.KbSpecificObject();
            String source = viccEntry.source();


            if (viccEntry.source().equals("brca")) {
                Brca kbBrca = (Brca) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("cgi")) {
                Cgi kbCgi = (Cgi) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("civic")) {
                Civic kbCivic = (Civic) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("jax")) {
                Jax kbJax = (Jax) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("jax_trials")) {
                JaxTrials kbJaxTrials = (JaxTrials) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("molecularmatch")) {
                MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("molecularmatch_trials")) {
                MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("oncokb")) {
                OncoKb kbOncoKb = (OncoKb) kbSpecificObject;

                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("pmkb")) {
                Pmkb kbPmkb = (Pmkb) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else if (viccEntry.source().equals("sage")) {
                Sage kbSage = (Sage) kbSpecificObject;
                String actionableCNV = "TODO";
                String knownCNV = "TODO";
            } else {
                LOGGER.warn("Unknown source");
            }
        }

    }
}
