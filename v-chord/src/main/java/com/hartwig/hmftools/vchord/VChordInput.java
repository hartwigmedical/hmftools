package com.hartwig.hmftools.vchord;

import ai.djl.modality.cv.Image;

public class VChordInput
{
    public final Image purpleCircos;
    public final HrdCancerType cancerType;
    public final double purity;

    public VChordInput(final Image purpleCircos, final HrdCancerType cancerType, final double purity)
    {
        this.purpleCircos = purpleCircos;
        this.cancerType = cancerType;
        this.purity = purity;
    }
}
