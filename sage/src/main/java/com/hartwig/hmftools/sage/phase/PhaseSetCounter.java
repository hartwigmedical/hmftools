package com.hartwig.hmftools.sage.phase;

public class PhaseSetCounter
{
    private int mPhaseSet;

    public PhaseSetCounter()
    {
        mPhaseSet = 1;
    }

    public synchronized int getNext()
    {
        return mPhaseSet++;
    }
}
