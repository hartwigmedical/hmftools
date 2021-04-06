package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.genome.position.GenomePositions.union;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PCFPositionsSupplier {

    private static final Logger LOGGER = LogManager.getLogger(PCFPositionsSupplier.class);

    @NotNull
    public static Multimap<Chromosome, PCFPosition> createPositions(@NotNull final ConfigSupplier configSupplier) {
        final CommonConfig config = configSupplier.commonConfig();
        final AmberData amberData = configSupplier.amberData();
        final CobaltData cobaltData = configSupplier.cobaltData();

        final Multimap<Chromosome, PCFPosition> referenceBreakPoint = cobaltData.referenceSegments();
        final Multimap<Chromosome, PCFPosition> tumorBreakPoints = cobaltData.tumorSegments();
        final Multimap<Chromosome, PCFPosition> tumorBAF = amberData.tumorSegments();

        LOGGER.info("Merging reference and tumor ratio break points");
        return union(union(tumorBreakPoints, referenceBreakPoint), tumorBAF);
    }

    //    @NotNull
    //    private static Multimap<Chromosome, PCFPosition> centromeres(int windowSize, @NotNull final Map<Chromosome, GenomePosition> centromeres) {
    //        final Multimap<Chromosome, PCFPosition> result = ArrayListMultimap.create();
    //        final Window window = new Window(windowSize);
    //        for (Map.Entry<Chromosome, GenomePosition> entry : centromeres.entrySet()) {
    //            final GenomePosition centromere = entry.getValue();
    //            result.put(entry.getKey(), create(centromere.chromosome(), window.start(centromere.position())));
    //        }
    //        return result;
    //    }

    @NotNull
    private static PCFPosition create(@NotNull final String chromosome, long position) {
        return ImmutablePCFPosition.builder()
                .chromosome(chromosome)
                .position(position)
                .source(PCFSource.REFERENCE_RATIO)
                .minPosition(position)
                .maxPosition(position)
                .build();
    }
}
