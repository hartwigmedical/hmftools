package com.hartwig.hmftools.bachelor.predicates;

import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.bachelor.VariantModel;

import nl.hartwigmedicalfoundation.bachelor.ProgramWhitelist;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WhitelistPredicate implements Predicate<VariantModel> {

    private final Multimap<String, String> transcriptProteinWhitelist = HashMultimap.create();
    private final Set<String> dbSNPWhitelist = Sets.newHashSet();
    private static final Logger LOGGER = LogManager.getLogger(WhitelistPredicate.class);

    public WhitelistPredicate(@NotNull final Multimap<String, String> geneToEnsemblMap, @Nullable final ProgramWhitelist whitelist) {
        if (whitelist != null) {
            for (final Object variantOrDbSNP : whitelist.getVariantOrDbSNP()) {
                if (variantOrDbSNP instanceof ProgramWhitelist.Variant) {
                    final ProgramWhitelist.Variant variant = (ProgramWhitelist.Variant) variantOrDbSNP;
                    for (final String transcript : geneToEnsemblMap.get(variant.getGene().getName())) {
                        transcriptProteinWhitelist.put(transcript, variant.getHGVSP());
                    }
                } else if (variantOrDbSNP instanceof String) {
                    dbSNPWhitelist.add((String) variantOrDbSNP);
                }
            }
        }
    }

    @Override
    public boolean test(final VariantModel variant) {
        if(inDbSNPWhitelist(variant))
        {
            LOGGER.debug("variant({}) found in dbSNP whitelist", variant.context().getID());
            return true;
        }
        else if(inProteinWhitelist(variant))
        {
            LOGGER.debug("variant({}) found in protein whitelist", variant.context().getID());
            return true;
        }

        return false;
    }

    private boolean inProteinWhitelist(final VariantModel variant) {
        return variant.sampleAnnotations()
                .stream()
                .filter(a -> !a.hgvsProtein().isEmpty())
                .filter(a -> transcriptProteinWhitelist.containsKey(a.featureID()))
                .anyMatch(a -> transcriptProteinWhitelist.get(a.featureID()).contains(a.hgvsProtein().replaceFirst("^p\\.", "")));
    }

    private boolean inDbSNPWhitelist(final VariantModel variant) {
        return variant.dbSNP().stream().anyMatch(dbSNPWhitelist::contains);
    }

}
