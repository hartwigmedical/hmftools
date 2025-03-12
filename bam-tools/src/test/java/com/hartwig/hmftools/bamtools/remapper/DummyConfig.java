package com.hartwig.hmftools.bamtools.remapper;

class DummyConfig extends AltContigRemapperConfig
{
    private final PairAligner mAligner;

    DummyConfig(PairAligner aligner, String inputBam, String outputBam, boolean overrideSliceHlaOnly)
    {
        super(inputBam, outputBam, null, overrideSliceHlaOnly);
        mAligner = aligner;
    }

    @Override
    public PairAligner pairAligner()
    {
        return mAligner;
    }
}
