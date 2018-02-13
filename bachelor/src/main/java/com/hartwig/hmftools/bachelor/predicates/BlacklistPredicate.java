package com.hartwig.hmftools.bachelor.predicates;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.SnpEff;
import com.hartwig.hmftools.bachelor.VariantModel;

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
        for (final SnpEff annotation : variantModel.annotations()) {
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

    private static boolean test(final ProgramBlacklist.Exclusion blacklist, final VariantContext context, final SnpEff annotation) {
        if (blacklist.getHGVSP() != null && !annotation.hgvsProtein().isEmpty() && blacklist.getHGVSP().equals(annotation.hgvsProtein())) {
            return true;
        }
        if (blacklist.getHGVSC() != null && !annotation.hgvsCoding().isEmpty() && blacklist.getHGVSC().equals(annotation.hgvsCoding())) {
            return true;
        }
        if (blacklist.getMinCodon() != null && !annotation.proteinPosition().isEmpty() // TODO: stronger check here?
                && blacklist.getMinCodon().intValue() <= annotation.proteinPosition().get(0)) {
            return true;
        }

        return blacklist.getPosition() != null && atPosition(context, blacklist.getPosition());
    }

    private static boolean atPosition(final VariantContext v, final String position) {
        // TODO: robust enough check?
        return position.equals(v.getContig() + ":" + v.getStart());

    }
}
