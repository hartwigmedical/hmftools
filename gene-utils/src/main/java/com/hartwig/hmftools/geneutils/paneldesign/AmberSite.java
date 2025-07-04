package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.BasePosition;

public class AmberSite extends BasePosition
{
    private ProbeCandidate mProbe;

    public AmberSite(final String chromosome, final int position)
    {
        super(chromosome, position);
        mProbe = null;
    }

    public ProbeCandidate probe()
    {
        return mProbe;
    }

    public void setProbe(final ProbeCandidate candidate)
    {
        mProbe = candidate;
    }
}
