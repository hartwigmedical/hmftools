package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

public class SomaticVariantData implements ComparableItem
{
    public final SomaticVariant Variant;

    public SomaticVariantData(final SomaticVariant variant)
    {
        Variant = variant;
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public boolean reportable()
    {
        return Variant.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SomaticVariantData otherVar = (SomaticVariantData)other;

        if(!Variant.chromosome().equals(otherVar.Variant.chromosome()) || Variant.position() != otherVar.Variant.position())
            return false;

        if(!Variant.ref().equals(otherVar.Variant.ref()) || !Variant.alt().equals(otherVar.Variant.alt()))
            return false;

        if(Variant.type() != otherVar.Variant.type())
            return false;

        return true;
    }

    @Override
    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final SomaticVariantData otherVar = (SomaticVariantData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Variant.reported(), otherVar.Variant.reported());

        if(matchLevel == REPORTABLE)
            return diffs;

        checkDiff(diffs, "hotspot", Variant.hotspot().toString(), otherVar.Variant.hotspot().toString());
        checkDiff(diffs, "canonicalEffect", Variant.canonicalEffect(), otherVar.Variant.canonicalEffect());
        checkDiff(diffs, "canonicalCodingEffect", Variant.canonicalCodingEffect().toString(), otherVar.Variant.canonicalCodingEffect().toString());
        checkDiff(diffs, "tier", Variant.tier().toString(), otherVar.Variant.tier().toString());
        checkDiff(diffs, "subclonalLikelihood", Variant.subclonalLikelihood(), otherVar.Variant.subclonalLikelihood());
        checkDiff(diffs, "variantCopyNumber", Variant.variantCopyNumber(), otherVar.Variant.variantCopyNumber());

        return diffs;
    }

    @Override
    public String description()
    {
        return String.format("%s:%d %s>%s %s",
                Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type());
    }
}
