package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_DISCORD_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_KNOWN_TYPE;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_REALIGN_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_SPLIT_FRAGS;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record RnaFusionData(RnaFusion RnaFusion, BasePosition ComparisonPositionUp,
                            BasePosition ComparisonPositionDown) implements ComparableItem
{
    public static String FLD_JUNC_TYPE_UP = "JuncTypeUp";
    public static String FLD_JUNC_TYPE_DOWN = "JuncTypeDown";

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%s:%d", RnaFusion.name(), RnaFusion.chromosomeUp(),
                RnaFusion.positionUp(), RnaFusion.chromosomeDown(), RnaFusion.positionDown());

        boolean upLifted = ComparisonPositionUp.Position != RnaFusion.positionUp()
                || !ComparisonPositionUp.Chromosome.equals(RnaFusion.chromosomeUp());

        boolean downLifted = ComparisonPositionDown.Position != RnaFusion.positionDown()
                || !ComparisonPositionDown.Chromosome.equals(RnaFusion.chromosomeDown());

        if(upLifted || downLifted)
            key += String.format(" liftover(%s-%s)", ComparisonPositionUp, ComparisonPositionDown);

        return key;
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", RnaFusion.knownType()));
        values.add(format("%s", RnaFusion.junctionTypeUp()));
        values.add(format("%s", RnaFusion.junctionTypeDown()));
        values.add(format("%d", RnaFusion.splitFragments()));
        values.add(format("%d", RnaFusion.realignedFrags()));
        values.add(format("%d", RnaFusion.discordantFrags()));

        return values;
    }

    @Override
    public boolean reportable()
    {
        return false;
    }

    @Override
    public boolean isPass()
    {
        return RnaFusion.filter().equals("PASS");
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaFusionData otherData = (RnaFusionData)other;

        if(!otherData.RnaFusion.name().equals(RnaFusion.name())){
            return false;
        }
        if(!otherData.RnaFusion.chromosomeUp().equals(ComparisonPositionUp.Chromosome))
        {
            return false;
        }
        if(!otherData.RnaFusion.chromosomeDown().equals(ComparisonPositionDown.Chromosome))
        {
            return false;
        }
        if(otherData.RnaFusion.positionUp() != ComparisonPositionUp.Position)
        {
            return false;
        }
        return otherData.RnaFusion.positionDown() == ComparisonPositionDown.Position;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final RnaFusion ref = RnaFusion;
        final RnaFusion otherData = ((RnaFusionData) other).RnaFusion;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_KNOWN_TYPE, ref.knownType().toString(), otherData.knownType().toString());
        checkDiff(diffs, FLD_JUNC_TYPE_UP, ref.junctionTypeUp(), otherData.junctionTypeUp());
        checkDiff(diffs, FLD_JUNC_TYPE_DOWN, ref.junctionTypeDown(), otherData.junctionTypeDown());
        checkDiff(diffs, FLD_SPLIT_FRAGS, ref.splitFragments(), otherData.splitFragments(), thresholds);
        checkDiff(diffs, FLD_REALIGN_FRAGS, ref.realignedFrags(), otherData.realignedFrags(), thresholds);
        checkDiff(diffs, FLD_DISCORD_FRAGS, ref.discordantFrags(), otherData.discordantFrags(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
