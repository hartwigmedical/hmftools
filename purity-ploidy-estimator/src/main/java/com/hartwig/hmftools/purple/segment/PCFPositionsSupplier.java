package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.position.GenomePositions.union;

import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.window.Window;
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
    public static Multimap<String, PCFPosition> createPositions(@NotNull final ConfigSupplier configSupplier) {
        final CommonConfig config = configSupplier.commonConfig();
        final AmberData amberData = configSupplier.amberData();
        final CobaltData cobaltData = configSupplier.cobaltData();

        final Multimap<String, PCFPosition> referenceBreakPoint = cobaltData.referenceSegments();
        final Multimap<String, PCFPosition> tumorBreakPoints = cobaltData.tumorSegments();
        final Multimap<String, PCFPosition> tumorBAF = amberData.tumorSegments();

        LOGGER.info("Merging reference and tumor ratio break points");
        return union(union(union(tumorBreakPoints, referenceBreakPoint), tumorBAF),
                centromeres(config.windowSize(), configSupplier.refGenomeConfig().centromere()));
    }

    @NotNull
    private static Multimap<String, PCFPosition> centromeres(int windowSize, @NotNull final Map<Chromosome, GenomePosition> centromeres) {
        final Multimap<String, PCFPosition> result = ArrayListMultimap.create();
        final Window window = new Window(windowSize);
        for (GenomePosition centromere : centromeres.values()) {
            result.put(centromere.chromosome(), create(centromere.chromosome(), window.start(centromere.position())));
        }

        return result;
    }

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
