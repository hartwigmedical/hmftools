package com.hartwig.hmftools.bachelor.predicates;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.bachelor.VariantModel;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;
import nl.hartwigmedicalfoundation.bachelor.ProgramWhitelist;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class WhitelistPredicate implements Predicate<VariantModel> {

    // private final Multimap<String, String> mTranscriptProteins = HashMultimap.create();

    private final ProgramWhitelist mWhitelist;

    private static final Logger LOGGER = LogManager.getLogger(WhitelistPredicate.class);

    public WhitelistPredicate(@NotNull final Multimap<String, String> geneToEnsemblMap, @Nullable final ProgramWhitelist whitelist)
    {
        mWhitelist = whitelist;

        /*
        if (whitelist != null)
        {
            for (final Object variantOrDbSNP : whitelist.getVariantOrDbSNP())
            {
                if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
                {
                    final ProgramWhitelist.Variant variant = (ProgramWhitelist.Variant) variantOrDbSNP;
                    for (final String transcript : geneToEnsemblMap.get(variant.getGene().getName()))
                    {
                        mTranscriptProteins.put(transcript, variant.getHGVSP());
                    }
                }
                else if (variantOrDbSNP instanceof String)
                {
                    mWhitelistDbSNPList.add((String) variantOrDbSNP);
                }
            }
        }
        */
    }

    @Override
    public boolean test(final VariantModel variantModel)
    {
        if(mWhitelist == null)
            return false;

        for (final SnpEffAnnotation annotation : variantModel.sampleAnnotations())
        {
            for (final Object variantOrDbSNP : mWhitelist.getVariantOrDbSNP())
            {
                if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
                {
                    final ProgramWhitelist.Variant whitelistVar = (ProgramWhitelist.Variant) variantOrDbSNP;

                    if(matchesWhitelistGeneProtein(whitelistVar, variantModel.context(), annotation))
                    {
                        return true;
                    }
                }
                else if (variantOrDbSNP instanceof String)
                {
                    if(matchesWhitelistDbSNPId((String) variantOrDbSNP, variantModel.context(), annotation))
                    {
                        return true;
                    }
                }
            }
        }

        return false;

        /*
            if (inDbSNPWhitelist(variant))
            {
                LOGGER.debug("variant({}) found in dbSNP whitelist", variant.context().getID());
                return true;
            }
            else if (inProteinWhitelist(variant))
            {
                LOGGER.debug("variant({}) found in protein whitelist", variant.context().getID());
                return true;
            }
     */
    }

    public static boolean matchesWhitelistGeneProtein(ProgramWhitelist.Variant geneProtein,
            final VariantContext context, final SnpEffAnnotation annotation)
    {
        if(!geneProtein.getGene().getName().equals(annotation.gene()))
            return false;

        if (geneProtein.getHGVSP() != null && !annotation.hgvsProtein().isEmpty()
        && geneProtein.getHGVSP().equals(annotation.hgvsProtein().replaceFirst("^p\\.", "")))
        {
            LOGGER.debug("variant({}) found in blacklist HGVSP({})", context.getID(), geneProtein.getHGVSP());
            return true;
        }

        return false;

        /*
        return variant.sampleAnnotations()
                .stream()
                .filter(a -> !a.hgvsProtein().isEmpty())
                .filter(a -> mTranscriptProteins.containsKey(a.transcript()))
                .anyMatch(a -> mTranscriptProteins.get(a.transcript()).contains(a.hgvsProtein().replaceFirst("^p\\.", "")));
        */
    }

    public static boolean matchesWhitelistDbSNPId(final String dnSNPId, final VariantContext variant, final SnpEffAnnotation annotation)
    {
        Set<String> varDbSNPList = Lists.newArrayList(variant.getID()
                .split(","))
                .stream().filter(s -> s.startsWith("rs"))
                .collect(Collectors.toSet());

        if(varDbSNPList.contains(dnSNPId))
        {
            LOGGER.debug("variant({}) matched in whitelist rs DB Ids list({})", variant.getID());
            return true;
        }
        else
        {
            return false;
        }
    }

    public static String asString(final Object variantOrDbSNP)
    {
        if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
        {
            final ProgramWhitelist.Variant geneProtein = (ProgramWhitelist.Variant) variantOrDbSNP;
            return String.format("gene(%s) protein(%s)", geneProtein.getGene().getName(), geneProtein.getHGVSP().toString());
        }
        else
        {
            return String.format("DbSNP ID(%s)", (String) variantOrDbSNP);
        }
    }

}
