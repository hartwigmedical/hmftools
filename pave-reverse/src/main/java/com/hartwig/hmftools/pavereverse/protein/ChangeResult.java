package com.hartwig.hmftools.pavereverse.protein;

import java.util.Objects;

import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;

public class ChangeResult
{
    public final AminoAcidSequence Acids;
    public final String ExonBases;
    final int Location;
    final String RefBases;
    final String AltBases;

    public ChangeResult(AminoAcidSequence aminoAcids, String exonBases, int location, String refBases, String altBases)
    {
        Acids = aminoAcids;
        ExonBases = exonBases;
        Location = location;
        RefBases = refBases;
        AltBases = altBases;
    }

    public BaseSequenceChange asChange(String chromosome)
    {
        return new BaseSequenceChange(RefBases, AltBases, chromosome, Location);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ChangeResult that = (ChangeResult) o;
        return Location == that.Location && Objects.equals(Acids, that.Acids)
                && Objects.equals(ExonBases, that.ExonBases) && Objects.equals(RefBases, that.RefBases)
                && Objects.equals(AltBases, that.AltBases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Acids, ExonBases, Location, RefBases, AltBases);
    }

    @Override
    public String toString()
    {
        return "ChangeResult{" +
                "AminoAcids=" + Acids +
                ", ExonBases='" + ExonBases + '\'' +
                ", Location=" + Location +
                ", RefBases='" + RefBases + '\'' +
                ", AltBases='" + AltBases + '\'' +
                '}';
    }
}
