package com.hartwig.hmftools.svtools.germline;

import java.util.List;

public class AlternatePath
{
    public final String VcfId;
    public final String MateId;
    public final List<Link> Path;

    public AlternatePath(final String vcfId, final String mateId, final List<Link> path)
    {
        VcfId = vcfId;
        MateId = mateId;
        Path = path;
    }
}
