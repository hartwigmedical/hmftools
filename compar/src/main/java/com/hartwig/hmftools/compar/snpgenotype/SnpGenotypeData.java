package com.hartwig.hmftools.compar.snpgenotype;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.compar.common.Category.SNP_GENOTYPE;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class SnpGenotypeData implements ComparableItem
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Genotype;
    public final BasePosition mComparisonPosition;

    protected static final String FLD_GENOTYPE = "Genotype";

    public SnpGenotypeData(final String chromosome, final int position, final String ref, final String alt, final String genotype,
            final BasePosition comparisonPosition)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Genotype = genotype;
        mComparisonPosition = comparisonPosition;
    }

    @Override
    public Category category()
    {
        return SNP_GENOTYPE;
    }

    @Override
    public String key()
    {
        if(mComparisonPosition.Position != Position)
        {
            return String.format("%s:%d %s liftover(%s)",
                    Chromosome, Position, Ref, mComparisonPosition);
        }
        else
        {
            return String.format("%s:%d %s", Chromosome, Position, Ref);
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", Alt));
        values.add(format("%s", Genotype));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SnpGenotypeData otherVar = (SnpGenotypeData) other;

        if(!mComparisonPosition.Chromosome.equals(otherVar.Chromosome) || mComparisonPosition.Position != otherVar.Position)
            return false;

        if(!Ref.equals(otherVar.Ref))
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final SnpGenotypeData otherData = (SnpGenotypeData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_ALT, Alt, otherData.Alt);
        checkDiff(diffs, FLD_GENOTYPE, Genotype, otherData.Genotype);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
