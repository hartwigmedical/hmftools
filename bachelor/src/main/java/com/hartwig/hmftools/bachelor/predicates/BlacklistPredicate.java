package com.hartwig.hmftools.bachelor.predicates;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.VariantModel;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class BlacklistPredicate implements Predicate<VariantModel> {

    @NotNull
    private final Collection<String> transcripts;
    @NotNull
    private final List<ProgramBlacklist.Exclusion> blacklist;

    private static final Logger LOGGER = LogManager.getLogger(BlacklistPredicate.class);

    public BlacklistPredicate(@NotNull final Collection<String> transcripts, @Nullable final ProgramBlacklist blacklist) {
        this.transcripts = transcripts;
        this.blacklist = blacklist != null ? blacklist.getExclusion() : Lists.newArrayList();
    }

    @Override
    public boolean test(final VariantModel variantModel) {
        for (final SnpEffAnnotation annotation : variantModel.sampleAnnotations()) {
            final boolean transcriptMatches = transcripts.contains(annotation.transcript());
            if (transcriptMatches) {
                for (ProgramBlacklist.Exclusion exclusion : blacklist) {
                    if (test(exclusion, variantModel.context(), annotation)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static boolean test(final ProgramBlacklist.Exclusion blacklist, final VariantContext context,
            final SnpEffAnnotation annotation) {

        if (blacklist.getHGVSP() != null && !annotation.hgvsProtein().isEmpty()
        && blacklist.getHGVSP().equals(annotation.hgvsProtein().replaceFirst("^p\\.", ""))) {

            LOGGER.debug("variant({}) found in blacklist HGVSP({})", context.getID(), blacklist.getHGVSP());
            return true;
        }
        if (blacklist.getHGVSC() != null && !annotation.hgvsCoding().isEmpty()
        && blacklist.getHGVSC().equals(annotation.hgvsCoding().replaceFirst("^c\\.", ""))) {

            LOGGER.debug("variant({}) found in blacklist HGVSC({})", context.getID(), blacklist.getHGVSC());

            return true;
        }

        final List<Integer> proteinPositions = proteinPosition(annotation);
        if (blacklist.getMinCodon() != null && !proteinPositions.isEmpty() // TODO: stronger check here?
        && blacklist.getMinCodon().intValue() <= proteinPositions.get(0)) {

            LOGGER.debug("variant({}) found in blacklist minCodon({})", context.getID(), blacklist.getMinCodon());
            return true;
        }

        if(blacklist.getPosition() != null && atPosition(context, blacklist.getPosition()))
        {
            LOGGER.debug("variant({}) found in blacklist postition({})", context.getID(), blacklist.getPosition());
            return true;
        }

        return false;
    }

    private static boolean atPosition(final VariantContext v, final String position) {
        // TODO: robust enough check?
        return position.equals(v.getContig() + ":" + v.getStart());
    }

    @NotNull
    private static List<Integer> proteinPosition(@NotNull final SnpEffAnnotation annotation) {
        return Arrays.stream(annotation.aaPosAndLength().split("/"))
                .filter(s -> !s.isEmpty())
                .map(Integer::parseInt)
                .collect(Collectors.toList());
    }
}
