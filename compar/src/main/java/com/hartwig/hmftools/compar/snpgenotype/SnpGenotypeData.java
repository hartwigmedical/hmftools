package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.compar.common.CategoryType.SNP_GENOTYPE;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class SnpGenotypeData implements ComparableItem
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Genotype;
    public final String VcfSampleId;
    public final BasePosition mComparisonPosition;

    public SnpGenotypeData(final String chromosome, final int position, final String ref, final String alt, final String genotype,
            final String vcfSampleId, final BasePosition comparisonPosition)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Genotype = genotype;
        VcfSampleId = vcfSampleId;
        mComparisonPosition = comparisonPosition;
    }

    @Override
    public CategoryType category()
    {
        return SNP_GENOTYPE;
    }

    @Override
    public String key()
    {
        if(mComparisonPosition.Position != Position)
        {
            return String.format("%s:%d %s liftover(%s)", Chromosome, Position, Ref, mComparisonPosition);
        }
        else
        {
            return String.format("%s:%d %s", Chromosome, Position, Ref);
        }
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
}
