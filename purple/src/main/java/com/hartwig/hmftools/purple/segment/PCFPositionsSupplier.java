package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.genome.position.GenomePositions.union;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PCFPositionsSupplier
{
    @NotNull
    public static Multimap<Chromosome, PCFPosition> createPositions(final AmberData amberData, final CobaltData cobaltData)
    {
        final Multimap<Chromosome, PCFPosition> referenceBreakPoint = cobaltData.ReferenceSegments;
        final Multimap<Chromosome, PCFPosition> tumorBreakPoints = cobaltData.TumorSegments;
        final Multimap<Chromosome, PCFPosition> tumorBAF = amberData.TumorSegments;

        PPL_LOGGER.info("Merging reference and tumor ratio break points");

        return union(union(tumorBreakPoints, referenceBreakPoint), tumorBAF);
    }
}
