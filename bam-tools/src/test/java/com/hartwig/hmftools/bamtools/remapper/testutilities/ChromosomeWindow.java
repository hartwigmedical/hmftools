package com.hartwig.hmftools.bamtools.remapper.testutilities;

import java.util.Arrays;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

record ChromosomeWindow(int chromosome, int start, int end)
{
    public ChromosomeWindow next(int step)
    {
        return new ChromosomeWindow(chromosome, start + step, end + step);
    }

    public String chromosomeName()
    {
        return "chr" + (chromosome + 1);
    }

    public BaseRegion toBaseRegion(RefGenomeSource refGenomeSource)
    {
        byte[] bases = refGenomeSource.getBases(chromosomeName(), start, end - 1);
        return new BaseRegion(chromosome, start, end, bases);
    }

    public BaseRegion mateBaseRegion(RefGenomeSource refGenomeSource)
    {
        int mateStart = end;
        int mateStop = mateStart + (end - start);
        byte[] mateBases = refGenomeSource.getBases(chromosomeName(), mateStart, mateStop - 1);
        return new BaseRegion(chromosome, mateStart, mateStop, mateBases);
    }
}
