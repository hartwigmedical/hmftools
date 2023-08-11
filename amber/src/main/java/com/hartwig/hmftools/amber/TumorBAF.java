package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class TumorBAF implements GenomePosition
{
    public int NormalReadDepth;
    public int NormalRefSupport;
    public int NormalAltSupport;

    public final PositionEvidence TumorEvidence;

    public TumorBAF(final String chromosome, final int position, final String ref, final String alt)
    {
        NormalReadDepth = 0;
        NormalRefSupport = 0;
        NormalAltSupport = 0;

        TumorEvidence = new PositionEvidence(chromosome, position, ref, alt);
    }

    @Override
    public String chromosome() { return TumorEvidence.Chromosome; }
    public int position() { return TumorEvidence.Position; }

    public String ref() { return TumorEvidence.ref(); }
    public String alt() { return TumorEvidence.alt(); }

    public double refFrequency() { return TumorEvidence.RefSupport / (double)TumorEvidence.ReadDepth; }
    public double altFrequency() {
        return TumorEvidence.AltSupport / (double)TumorEvidence.ReadDepth;
    }

    public static TumorBAF fromNormal(final PositionEvidence normal)
    {
        TumorBAF tumorBAF = new TumorBAF(normal.Chromosome, normal.Position, normal.ref(), normal.alt());
        tumorBAF.NormalReadDepth = normal.ReadDepth;
        tumorBAF.NormalRefSupport = normal.RefSupport;
        tumorBAF.NormalAltSupport = normal.AltSupport;
        return tumorBAF;
    }
}
