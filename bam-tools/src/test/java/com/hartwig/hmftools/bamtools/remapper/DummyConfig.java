package com.hartwig.hmftools.bamtools.remapper;

class DummyConfig extends AltContigRemapperConfig
{
    private final PairAligner mAligner;

    DummyConfig(PairAligner aligner, String inputBam, String outputBam)
    {
        super(inputBam, outputBam, null);
        mAligner = aligner;
    }

    @Override
    public PairAligner aligner()
    {
        return mAligner;
    }
}
