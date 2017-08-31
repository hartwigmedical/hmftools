package com.hartwig.hmftools.breakpointinspector.clipping;

import com.hartwig.hmftools.breakpointinspector.Location;

public class ClipInfo {
    Location Alignment;
    int Length = 0;
    String Sequence = "";
    boolean HardClipped = false;
    boolean Left = false;
    boolean Right = false;
}
