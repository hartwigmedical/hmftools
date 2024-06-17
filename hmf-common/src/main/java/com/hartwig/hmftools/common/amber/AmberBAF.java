package com.hartwig.hmftools.common.amber;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.region.BasePosition;

public class AmberBAF implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final double TumorBAF;
    public final int TumorDepth;
    public final double NormalBAF;
    public final int NormalDepth;

    public AmberBAF(
            final String chromosome, final int position, final double tumorBAF, final int tumorDepth, final double normalBAF,
            final int normalDepth)
    {
        Chromosome = chromosome;
        Position = position;
        TumorBAF = tumorBAF;
        TumorDepth = tumorDepth;
        NormalBAF = normalBAF;
        NormalDepth = normalDepth;
    }

    @Override
    public String chromosome() { return Chromosome; }
    public int position() { return Position; }

    public double tumorModifiedBAF()
    {
        return 0.5 + Math.abs(tumorBAF() - 0.5);
    }
    public double normalModifiedBAF()
    {
        return 0.5 + Math.abs(normalBAF() - 0.5);
    }

    public String toString()
    {
        return format("location(%s) tumor(baf=%.4f depth=%d) normal(baf=%.4f depth=%d)",
                    super.toString(), TumorBAF, TumorDepth, NormalBAF, NormalDepth);
    }

    // backwards compatibility
    public double tumorBAF() { return TumorBAF; }
    public int tumorDepth() { return TumorDepth; }
    public double normalBAF() { return NormalBAF; }
    public int normalDepth() { return NormalDepth; }
}
