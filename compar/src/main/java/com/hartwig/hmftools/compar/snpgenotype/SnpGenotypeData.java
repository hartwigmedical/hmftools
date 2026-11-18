package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.compar.common.CategoryType.SNP_GENOTYPE;

import java.util.List;

import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class SnpGenotypeData extends ComparableItem
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Genotype;
    public final String VcfSampleId;

    public SnpGenotypeData(
            final String chromosome, final int position, final String ref, final String alt, final String genotype,
            final String vcfSampleId, final List<FieldInfo> fields)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Genotype = genotype;
        VcfSampleId = vcfSampleId;

        addStringValue(SnpGenotypeComparer.Fields.Alt.toString(), alt, fields);
        addStringValue(SnpGenotypeComparer.Fields.Genotype.toString(), genotype, fields);
        addStringValue(SnpGenotypeComparer.Fields.VcfSampleId.toString(), vcfSampleId, fields);
    }

    @Override
    public CategoryType category()
    {
        return SNP_GENOTYPE;
    }

    @Override
    public String key()
    {
        return String.format("%s:%d %s", Chromosome, Position, Ref);
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SnpGenotypeData otherVar = (SnpGenotypeData) other;

        if(!Chromosome.equals(otherVar.Chromosome) || Position != otherVar.Position)
            return false;

        if(!Ref.equals(otherVar.Ref))
            return false;

        return true;
    }
}
