package com.hartwig.hmftools.bachelor.predicates;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.VariantModel;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class BlacklistPredicate implements Predicate<VariantModel> {

    @NotNull
    private final Collection<String> transcripts;
    @NotNull
    private final List<ProgramBlacklist.Exclusion> blacklist;

    public BlacklistPredicate(@NotNull final Collection<String> transcripts, @Nullable final ProgramBlacklist blacklist) {
        this.transcripts = transcripts;
        this.blacklist = blacklist != null ? blacklist.getExclusion() : Lists.newArrayList();
    }

    @Override
    public boolean test(final VariantModel variantModel) {
        for (final VariantAnnotation annotation : variantModel.sampleAnnotations()) {
            final boolean transcriptMatches = transcripts.contains(annotation.featureID());
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
            final VariantAnnotation annotation) {
        if (blacklist.getHGVSP() != null && !annotation.hgvsProtein().isEmpty() && blacklist.getHGVSP()
                .equals(annotation.hgvsProtein().replaceFirst("^p\\.", ""))) {
            return true;
        }
        if (blacklist.getHGVSC() != null && !annotation.hgvsCoding().isEmpty() && blacklist.getHGVSC()
                .equals(annotation.hgvsCoding().replaceFirst("^c\\.", ""))) {
            return true;
        }

        final List<Integer> proteinPositions = proteinPosition(annotation);
        if (blacklist.getMinCodon() != null && !proteinPositions.isEmpty() // TODO: stronger check here?
                && blacklist.getMinCodon().intValue() <= proteinPositions.get(0)) {
            return true;
        }

        return blacklist.getPosition() != null && atPosition(context, blacklist.getPosition());
    }

    private static boolean atPosition(final VariantContext v, final String position) {
        // TODO: robust enough check?
        return position.equals(v.getContig() + ":" + v.getStart());
    }

    @NotNull
    private static List<Integer> proteinPosition(@NotNull final VariantAnnotation annotation) {
        return Arrays.stream(annotation.aaPosAndLength().split("/"))
                .filter(s -> !s.isEmpty())
                .map(Integer::parseInt)
                .collect(Collectors.toList());
    }
}
