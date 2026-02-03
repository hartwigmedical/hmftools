package com.hartwig.hmftools.common.purple;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GermlineDeletionTestFactory
{
    public static GermlineAmpDel create(final String geneName)
    {
        return create(geneName, false, GermlineStatus.HET_DELETION, 0D, 0, 0);
    }

    public static GermlineAmpDel create(final String geneName, boolean reported, final String chromosome,
            final String chromosomeBand)
    {
        return create(geneName, reported, GermlineStatus.HET_DELETION, 0D, 0, 0, chromosome, chromosomeBand);
    }

    public static GermlineAmpDel create(final String geneName, boolean reported, final GermlineStatus tumorStatus)
    {
        return create(geneName, reported, tumorStatus, 0D, 0, 0);
    }

    public static GermlineAmpDel create(
            final String geneName, boolean reported, final GermlineStatus tumorStatus, double tumorCopyNumber)
    {
        return create(geneName, reported, tumorStatus, tumorCopyNumber, 0, 0);
    }

    public static GermlineAmpDel create(
            final String geneName, boolean reported, final GermlineStatus tumorStatus, double tumorCopyNumber, int regionStart, int regionEnd)
    {
        return create(geneName, reported, tumorStatus, tumorCopyNumber, regionStart, regionEnd, Strings.EMPTY, Strings.EMPTY);
    }

    public static GermlineAmpDel create(
            final String geneName, boolean reported, final GermlineStatus tumorStatus,
            double tumorCopyNumber, int regionStart, int regionEnd, final String chromosome, final String chromosomeBand)
    {
        return new GermlineAmpDel(
                geneName, "",
                chromosome,
                chromosomeBand,
                regionStart,
                regionEnd,
                0,
                0,
                0,
                GermlineDetectionMethod.SEGMENT,
                tumorStatus,
                tumorStatus,
                0D,
                tumorCopyNumber,
                Strings.EMPTY,
                0,
                reported ? ReportedStatus.REPORTED : ReportedStatus.NONE);
    }
}
